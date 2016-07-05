#ifndef _INITIALS_H
#define _INITIALS_H

void print_greetings(void);
void print_usage(int status);
void print_version(int status);
void generate_template_config(int status);
void parse_config(FILE *file);
 int parse_channels(FILE *file, CHA ch_tab[]);
void parse_options(int argc, char *argv[], char **f);

#endif
