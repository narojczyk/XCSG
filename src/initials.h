#ifndef _INITIALS_H
#define _INITIALS_H

void print_greetings(void);
void print_help(int status);
void print_info(int status);
void print_version(int status);
void generate_template_config(int status);

 int parse_config(FILE *file, CONFIG *cfg);
 int parse_config_legacy(FILE *file, CONFIG *cfg);
 int parse_channels(FILE *file, CHA ch_tab[], int num_channels);
 int parse_slits(FILE *file, SLI sl_tab[], int num_slits);
void parse_options(int argc, char *argv[], char **f);

#endif
