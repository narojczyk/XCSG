#ifndef _INITIALS_H
#define _INITIALS_H

void print_greetings(void);
void print_help(int status);
void print_info(int status);
void print_version(int status);
void generate_template_config(int status);

 int parse_config(FILE *file, CONFIG *cfg);
 int parse_config_legacy(FILE *file, CONFIG *cfg);
INC* parse_inclusions(FILE *file);
void parse_options(int argc, char *argv[], char **f);

#endif
