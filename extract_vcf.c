#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

// Function to print usage information
void print_usage(char *program_name) {
    fprintf(stderr, "Usage: %s [--id] [--info <INFO_FIELDS>] [--format <FORMAT_FIELDS>] [--split-fields <FIELDS_TO_SPLIT>] [--delimiter <DELIMITER>] <input.vcf> <output.tsv>\n", program_name);
    fprintf(stderr, "Example: %s --id --info AC,AF --format GT,DP input.vcf output.tsv\n", program_name);
}

// Function to parse comma-separated fields into an array
int parse_fields(char *fields_str, char **fields_arr[]) {
    char **fields = malloc(0);
    int count = 0;
    char *token = strtok(fields_str, ",");
    while (token != NULL) {
        count++;
        fields = realloc(fields, count * sizeof(char *));
        fields[count-1] = strdup(token);
        token = strtok(NULL, ",");
    }

    *fields_arr = fields;
    return count;
}

void split_fields(char delimiter, const char *fields, char ***field_array, int *size) {
    char *dstr = strdup(fields);

    // Calculate the number of substrings and allocate
    int count = 1;
    for (const char *tmp = dstr; *tmp != '\0'; tmp++) {
        if (*tmp == delimiter) {
            count++;
        }
    }
    *field_array = malloc(count * sizeof(char *));

    // Split string
    int i = 0;
    char *start = dstr;
    char *end = strchr(start, delimiter);
    if (end == NULL) { // Handle case where string cannot be split
        (*field_array)[i++] = strdup(start);
    }

    while (end != NULL) {
        size_t length = end - start;
        if (length == 0) continue; // String has trailing delimiter

        (*field_array)[i] = malloc((length) * sizeof(char));
        strncpy((*field_array)[i], start, length);
        (*field_array)[i][length] = '\0';
        i++;
        start = end + 1;
        end = strchr(start, delimiter);
    }

    *size = i;

    free(dstr);

}

void free_split_fields(char **field_array, int size) {
    for (int i = 0; i < size; i++) {
        free(field_array[i]);
    }
    free(field_array);
}

// Function to check if a field name exists in header
int field_exists_in_header(bcf_hdr_t *hdr, char *field_name) {
    return bcf_hdr_id2int(hdr, BCF_DT_ID, field_name) >= 0;
}

void put_info_value(bcf_hdr_t *hdr, bcf1_t *rec, char *tag, kstring_t *s) {
    kputc('\t', s);

    if (rec->d.info == NULL) {
        kputc('.', s);
        return;
    }

    bcf_info_t *info = bcf_get_info(hdr, rec, tag);
    if (info == NULL) {
        kputc('.', s);
        return;
    }

    if (info->len != 1) {
        bcf_fmt_array(s, info->len, info->type, info->vptr);
        return;
    }

    if (info->type == BCF_BT_FLOAT) {
        if (bcf_float_is_missing(info->v1.f) )
            kputc_('.', s);
        else
            kputd(info->v1.f, s);
    } else if (info->type == BCF_BT_CHAR) {
        kputc_(info->v1.i, s);
    } else if (info->type <= BCF_BT_INT32) {
        int64_t missing[] = {
            0, // BCF_BT_NULL
            bcf_int8_missing,
            bcf_int16_missing,
            bcf_int32_missing,
        };
        if (info->v1.i == missing[info->type])
            kputc_('.', s);
        else
            kputw(info->v1.i, s);
    }

    return;
}

void put_format_value(bcf_hdr_t *hdr, bcf1_t *rec, char *tag, int sample, kstring_t *s) {
    kputc('\t', s);

    if (rec->d.fmt == NULL) {
        kputc('.', s);
        return;
    }

    bcf_fmt_t *fmt = bcf_get_fmt(hdr, rec, tag);
    if (fmt == NULL) {
        kputc('.', s);
        return;
    }

    if (strcmp(tag, "GT") == 0) {
        bcf_format_gt(fmt, sample, s);
    } else {
        bcf_fmt_array(s, fmt->n, fmt->type, fmt->p + sample * (size_t) fmt->size);
    }
}

// Function to get INFO field value
void get_info_value(bcf_hdr_t *hdr, bcf1_t *rec, int info_idx, char *field_name, FILE *out_fp) {
    if (rec->d.info == NULL) {
        fprintf(out_fp, ".");
        return;
    }

    bcf_info_t *info = bcf_get_info(hdr, rec, field_name);
    if (info == NULL) {
        fprintf(out_fp, ".");
        return;
    }

    if (info->type == BCF_BT_CHAR && info->len == 1) {
        fprintf(out_fp, "%c", info->v1.i);
    } else if (info->type == BCF_BT_CHAR && info->len > 1) {
        fprintf(out_fp, "%.*s", info->len, info->vptr);
    } else if (info->type == BCF_BT_FLOAT) {
        fprintf(out_fp, "%.6f", info->v1.f);
    } else if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32) {
        fprintf(out_fp, "%d", info->v1.i);
    } else {
        fprintf(out_fp, ".");
    }
}

// Function to get FORMAT field value
void get_format_value(bcf_hdr_t *hdr, bcf1_t *rec, char *field_name, FILE *out_fp) {
    if (rec->d.fmt == NULL) {
        fprintf(out_fp, ".");
        return;
    }

    int format_idx = -1;
    for(int i = 0; i < rec->n_fmt; i++) {
        if (strcmp(hdr->id[BCF_DT_ID][rec->d.fmt[i].id].key, field_name) == 0) {
            format_idx = i;
            break;
        }
    }

    if (format_idx == -1) {
        fprintf(out_fp, ".");
        return;
    }

    bcf_fmt_t *fmt = bcf_get_fmt(hdr, rec, field_name);
    for (int j = 0; j < fmt->n; j++) {
        if (j > 0) {
            fprintf(out_fp, ",");
        }
        if (fmt->type == BCF_BT_CHAR && fmt->size == 1) {
            fprintf(out_fp, "%c", fmt->p[j]);
        } else if (fmt->type == BCF_BT_CHAR && fmt->size > 1) {
            fprintf(out_fp, "%.*s", fmt->size, &fmt->p[j * fmt->size]);
        } else if (fmt->type == BCF_BT_FLOAT) {
            float value;
            memcpy(&value, &fmt->p[j * sizeof(float)], sizeof(float));
            fprintf(out_fp, "%.6f", value);
        } else if (fmt->type == BCF_BT_INT8 || fmt->type == BCF_BT_INT16 || fmt->type == BCF_BT_INT32) {
            int32_t value;
            memcpy(&value, &fmt->p[j * sizeof(int32_t)], sizeof(int32_t));
            fprintf(out_fp, "%d", value);
        } else {
            fprintf(out_fp, ".");
        }
    }
}

int main(int argc, char *argv[]) {
    // Check the number of arguments
    if (argc < 4) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    // Parse command line arguments
    char *vcf_file = NULL;
    char *output_file = NULL;
    int include_id = 0;
    char *info_fields_str = NULL;
    char *format_fields_str = NULL;
    char *split_fields_str = NULL;
    char delimiter = ',';

    for (int i = 1; i < argc - 2; i++) {
        if (strcmp(argv[i], "--id") == 0) {
            include_id = 1;
        } else if (strcmp(argv[i], "--info") == 0 && i + 1 < argc - 2) {
            info_fields_str = argv[++i];
        } else if (strcmp(argv[i], "--format") == 0 && i + 1 < argc - 2) {
            format_fields_str = argv[++i];
        } else if (strcmp(argv[i], "--split-fields") == 0 && i + 1 < argc - 2) {
            split_fields_str = argv[++i];
        } else if(strcmp(argv[i], "--delimiter") == 0 && i + 1 < argc - 2) {
            delimiter = *(argv[++i]);
        } else {
            fprintf(stderr, "Error: Unknown option or missing argument: %s\n", argv[i]);
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    vcf_file = argv[argc - 2];
    output_file = argv[argc - 1];

    // Open VCF file for reading
    htsFile *fp = hts_open(vcf_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: Failed to open input VCF file: %s\n", vcf_file);
        return EXIT_FAILURE;
    }

    // Open output TSV file for writing
    FILE *out_fp = fopen(output_file, "w");
    if (out_fp == NULL) {
        fprintf(stderr, "Error: Failed to open output TSV file for writing: %s\n", output_file);
        hts_close(fp);
        return EXIT_FAILURE;
    }

    // Read VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (hdr == NULL) {
        fprintf(stderr, "Error: Failed to read VCF header\n");
        fclose(out_fp);
        hts_close(fp);
        return EXIT_FAILURE;
    }

    // Parse INFO and FORMAT fields if provided
    char **info_fields = NULL;
    int num_info_fields = 0;
    if (info_fields_str != NULL) {
        num_info_fields = parse_fields(info_fields_str, &info_fields);
    }

    char **format_fields = NULL;
    int num_format_fields = 0;
    if (format_fields_str != NULL) {
        num_format_fields = parse_fields(format_fields_str, &format_fields);
    }

    // Parse subset of INFO fields to split
    char **split_info_fields = NULL;
    int num_split_fields = 0;
    if (split_fields_str != NULL) {
        num_split_fields = parse_fields(split_fields_str, &split_info_fields);
        for (int i = 0; i < num_info_fields; i++) {
            for (int j = 0; j < num_split_fields; j++) {
                if (strcmp(info_fields[i], split_info_fields[j]) == 0 ) {
                    continue;
                }
            } 
        }

    }

    // Write header line to output file
    fprintf(out_fp, "CHROM\tPOS\tREF\tALT");
    if (include_id) {
        fprintf(out_fp, "\tID");
    }
    for (int i = 0; i < num_info_fields; i++) {
        fprintf(out_fp, "\t%s", info_fields[i]);
    }
    for (int i = 0; i < num_format_fields; i++) {
        fprintf(out_fp, "\t%s", format_fields[i]);
    }
    
    // Check if a FORMAT field exists
    int nsamples = bcf_hdr_nsamples(hdr);
    if (nsamples > 0 && num_format_fields > 0) fprintf(out_fp, "\tSAMPLE");
    fprintf(out_fp, "\n");

    // Iterate through variants and write to output file
    bcf1_t *rec = bcf_init();
    while (bcf_read(fp, hdr, rec) == 0) {
        kstring_t s = {0, 0, 0};

        bcf_unpack(rec, BCF_UN_ALL);

        kputs(bcf_hdr_id2name(hdr, rec->rid), &s); // CHROM
        kputc_('\t', &s); kputl(rec->pos + 1, &s); // POS
        kputc_('\t', &s); kputs(rec->d.allele[0], &s); // REF
        kputc_('\t', &s); // ALT

        if (rec->n_allele > 1) {
            for (int i = 1; i < rec->n_allele; ++i) {
                if (i > 1) kputc_(',', &s);
                kputs(rec->d.allele[i], &s);
            }
        } else kputc_('.', &s);

        if (include_id) { kputc_('\t', &s); kputs(rec->d.id, &s); } // ID

        for (int i = 0; i < num_info_fields; i++) put_info_value(hdr, rec, info_fields[i], &s);

        for (int i = 0; i < num_split_fields; i++) {
            char **fields_array = NULL;
            int num_fields = 0;
            split_fields(delimiter, split_info_fields[i], &fields_array, &num_fields);

            for (int j = 0; j < num_fields; j++) printf("%s\t", fields_array[j]);
        }
        printf("\n");

        if (nsamples == 0 || num_format_fields == 0) {
            fprintf(out_fp, "%s\n", s.s);
            continue;
        }

        for (int n = 0; n < nsamples; n++) { // If format field exists
            kstring_t ss = {0};
            kputs(s.s, &ss);
            for (int i = 0; i < num_format_fields; i++) put_format_value(hdr, rec, format_fields[i], n, &ss);

            kputc('\t', &ss); kputs(hdr->samples[n], &ss); // SAMPLE NAME
            fprintf(out_fp, "%s\n", ss.s);

            free(ss.s);
        }
        free(s.s);
    }

    // Clean up
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    fclose(out_fp);

    // Free allocated memory for fields
    for (int i = 0; i < num_info_fields; i++) {
        free(info_fields[i]);
    }
    free(info_fields);

    for (int i = 0; i < num_format_fields; i++) {
        free(format_fields[i]);
    }
    free(format_fields);

    return EXIT_SUCCESS;
}
