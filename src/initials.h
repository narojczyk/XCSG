#ifndef _INITIALS_H
#define _INITIALS_H

void print_greetings(void);
void print_help(int status);
void print_info(int status);
void print_version(int status);
void generate_template_config(int status);
 int parse_config(FILE *file);
 int parse_channels(FILE *file, CHA ch_tab[]);
 int parse_slits(FILE *file, SLI sl_tab[]);
void parse_options(int argc, char *argv[], char **f);

#endif
