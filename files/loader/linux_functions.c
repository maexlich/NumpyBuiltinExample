#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include <errno.h>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include "linker_generator.h"
//#include <glib.h>

int getExecPath_l(char *path, size_t len, int get_dirname){
	if (readlink("/proc/self/exe", path, len) != -1) {
        if(get_dirname)
            dirname(path); 
        return TRUE;
    }
    else
    	return FALSE;
}

int createFailIfExists_l(char *file){
	int fd;
	fd = open(file, O_CREAT|O_EXCL|O_RDWR, S_IRUSR|S_IWUSR);
    if (fd == -1 || errno == EEXIST)
    	return -1;	// File exists
    else {
    	close(fd);
    	return 0;	// File created
    }
}

int DeleteDirectory_l(const char *path, int nothing)
{
   DIR *d = opendir(path);
   size_t path_len = strlen(path);
   int r = -1;

   if (d)
   {
      struct dirent *p;

      r = 0;

      while (!r && (p=readdir(d)))
      {
          int r2 = -1;
          char *buf;
          size_t len;

          /* Skip the names "." and ".." as we don't want to recurse on them. */
          if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
          {
             continue;
          }

          len = path_len + strlen(p->d_name) + 2; 
          buf = malloc(len);

          if (buf)
          {
             struct stat statbuf;

             snprintf(buf, len, "%s/%s", path, p->d_name);

             if (!stat(buf, &statbuf))
             {
                if (S_ISDIR(statbuf.st_mode))
                {
                   r2 = DeleteDirectory_l(buf, 1);
                }
                else
                {
                   r2 = unlink(buf);
                }
             }

             free(buf);
          }

          r = r2;
      }

      closedir(d);
   }

   if (!r)
   {
      r = rmdir(path);
   }

   return r;
}