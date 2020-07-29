    #include <dirent.h>  
    #include <stdio.h> 
    #include <string.h> 
    #include <stdlib.h>
    int main(void)
    { FILE *fp;
      DIR           *d;
      struct dirent *dir;
      char *filenames;
      int a;
      int i = 0;
      d = opendir("revisions");
      if (d)
      {
        while ((dir = readdir(d)) != NULL)
        {
          filenames[i] = malloc(strlen(dir->d_name)+1);
          strcpy(filenames[i],dir->d_name);
          printf("%s", filenames);
          printf("\n");
          i++;
        }
        printf("\n");
        closedir(d);
      }
      for (a = 2; a < 8; a++){
         printf("%s", filenames[a]);
         printf("\n");
      }
      return(0);
    }
