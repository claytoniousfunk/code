
 #include <dirent.h>  
 #include <stdio.h> 
 #include <string.h> 
 #include <stdlib.h>
  
void openDir() 
{ 
    struct dirent *de;  // Pointer for directory entry 
    
  
    // opendir() returns a pointer of DIR type.  
    DIR *dr = opendir("/home/clayton/Analysis/skims"); 

   if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    { 
        printf("Could not open current directory" ); 
        return 0; 
    } 
	// count files in directory
	/*int file_count = 0;
	struct dirent * entry;
	while ((entry = readdir(dr)) != NULL) {
    	if (entry->d_type == DT_REG) { 
         	file_count++;
    	}
	}
	closedir(dr);
  
  printf("%i\n", file_count);
  */
  char filenames[256][256]; 
   
   
    int i = 0;
    if(dr){   
	while((de=readdir(dr)) != NULL) {
            //printf("%s\n", de->d_name); 
	  //filenames[i] = static_cast<char*>(malloc(strlen(de->d_name[i])+1));
	  
	  //filenames[i] = (char*) malloc(sizeof(de->d_name));
          strcpy(filenames[i],de->d_name);
          //printf("%s", filenames[i]);
          //printf("\n");
	  i++;
    	    }
 	closedir(dr);
	}

printf("%s\n", filenames[14]);



    return; 
} 

