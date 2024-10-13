//Confero is a utility that uses statistics to search a set of files for ones that are similar. Like different edits of the same document, or different archives which contain the same files.
//Version 1.0b. First public release. May contain bugs.
//Released under the GNU General Public Licence V3.


//Due to use of memory-mapped files, compiling for Windows needs an extra .c file and alternate mman.h that wraps the Windows memory map functions.
// gcc confero.c -O3 -o confero -Wall
// i686-w64-mingw32-gcc confero.c  -O3 -o confero.exe -Wall mman.c


#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <sys/stat.h>
#ifdef __MINGW32__
  #include "mman.h" //Alternate version that wraps Windows memory map functions.
#else
  #include <sys/mman.h>
#endif
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <inttypes.h>


int hash(char *filename, uint8_t *bloomhash);
uint32_t get_chunk_len(uint8_t *data, uint32_t max, uint8_t threshhold);
uint64_t make_fnv_from_block(uint8_t *block, uint32_t blocksize);
uint32_t bloom_bytes = 1024*8; //The actual number of buckets is eight times this. Must be a multiple of 8 (due to __builtin_popcountll)
void diag_show_hash(uint8_t *hash);
float jaccard_simularity(uint64_t *a, uint64_t *b);
void dofile(char *filename, struct stat *statbuf);
void dofolder( char *foldername);
uint64_t get_file_size(char *filename);
uint8_t *bloomhashes=NULL;
char **filenames=NULL;
int num_files=0;
uint64_t alloc_files=0;
float match_thresh=0.9;
uint8_t chunking_threshhold=52; //Lower number, bigger chunks. 52 seems about right.
//#define diag_mode 1
#define alloc_chunk 256

int main(int argc, char *argv[]){
  uint32_t n;
  for(n=1;n<argc;n++){
    if(argv[n][0] == '-'){
      if(strcmp(argv[n],"--help") == 0){
        printf("Confero: A file comparison program. Through the use of variable length chunking, Bloom filters and Jaccard similarity metric, compares a set of files to find similar ones.\n");
        printf("         Specifically it looks for long strings of data in common between two or more files within a (potentially very large) set.\n");
        printf("         For example, it will find different edits of a document that all descend from a common source. Or versions with different metadate.\n");
        printf("         Most usefully, it does this regardless of format. Text, images, audio, it matters not - it's all just data. Format-agnostic. Though compression may throw it off.\n");
        printf("         Only a problem if two versions use different types of settings of compression, as it'lll be working on the compressed byte stream.\n\n");
        printf("Usage:   confero [options] <file> [file] [file] ...\n");
        printf("         If you specify a folder it'll process that folder recursively.\n\n");
        printf("Options: -Bn            Set number of buckets. Default is -B1024. More buckets means more accurate comparisons, at the expense of needing more RAM.\n");
        printf("                        Needs buckets * 8 bytes per file. So default needs 8KiB per file.\n");
        printf("         -Tn            Chunking threshhold. Default -T52. Higher value, more chunks, will identify smaller matches and so more accurate comparisons. But also need more buckets.\n");
        printf("                        Sensible values are 47-55. Note that too small a value won't produce enough data to match on small files - you'll get a warning if this happens.\n");
        printf("                        If you receive warnings about files being too large to process, increase one or both of these B or T.\n");
        printf("         -Mn            Percentage of simularity to consider a match. Any file-pair with a simularity equal ot greater than this will be output.\n\n");
        printf("         Confero is a *statistical* utility. It does not produce guaranteed results. False positives are possible.\n");
        printf("         Setting B too high will never be harmful, except to memory usage. Increasing it is required to process very large files, and increases accuracy at the expense of memory.\n");
        printf("         Too low a T will fail on small files, too high will fail on large files unless B is increased to compensate.\n");
        printf("         Should either of these situations occur, a warning will be output.\n");
        return(1);
      }
     int mt;
      switch(argv[n][1]){
        case 'B':
          bloom_bytes=atoi(argv[n]+2)*8;
          if(!bloom_bytes){
             printf("Invalid B value.\n");
             return(1);
           }
          break;
        case 'T':
          chunking_threshhold=atoi(argv[n]+2);
          if(!chunking_threshhold){
            printf("Invalid T value.\n");
            return(1);
          }
          if(chunking_threshhold>60 ||chunking_threshhold < 45){
            printf("Sensible values of T are 45-60. This is an exponenent, so a single step in value will double or halve average chunk size.\n");
            return(1);
          }
          break;
        case 'M':
          mt=atoi(argv[n]+2);
          if(mt<1 || mt > 100 ){
             printf("Invalid M value.\n");
             return(1);
          }
          match_thresh=(float)mt / 100;
          break;
        default:
          printf("Unknown option %s\n", argv[n]);
          return(1);
      }
    }
  }


  for(n=1;n<argc;n++){
    if(argv[n][0]!= '-'){
      int arglen=strlen(argv[n]);
      if(argv[n][arglen-1] == '/')
         argv[n][arglen-1] = 0;
      dofolder(argv[n]);
    }
  }
  printf("Read and processed %u files.\nB=%u T=%u\n", num_files, bloom_bytes/8, chunking_threshhold);
  uint32_t a,b;
  #ifdef diag_mode
    for(a=0;a<num_files;a++){
      printf("%s\n", filenames[a]);
      diag_show_hash(bloomhashes+(a*bloom_bytes));
    }
  #endif
  for(a=0;a<num_files;a++){
    for(b=0;b<a;b++){
      float simularity=jaccard_simularity((uint64_t*)(bloomhashes+(a*bloom_bytes)),(uint64_t*)(bloomhashes+(b*bloom_bytes)));
      if(simularity >= match_thresh){
        printf("%f\n  %s\n  %s\n", simularity, filenames[a], filenames[b]);
        printf("%" PRIu64  ",%" PRIu64 "\n\n", get_file_size(filenames[a]),get_file_size(filenames[b]));
      }
    }
  }
}

float jaccard_simularity(uint64_t *a, uint64_t *b){ //Will need to cast the pointer when calling this.
  uint32_t comparisons=bloom_bytes/8;
  uint32_t n;
  uint32_t jac_union=0;
  uint32_t jac_intersection=0;
  for(n=0;n<comparisons;n++){
    jac_union+=__builtin_popcountll(a[n] | b[n]);
    jac_intersection+=__builtin_popcountll(a[n] & b[n]);
  }
  if(jac_union==0 || jac_intersection ==0)
    return(-1);
  return((float)jac_intersection / (float)jac_union);
}

void diag_show_hash(uint8_t *hash){
  uint32_t n;
  for(n=0;n<bloom_bytes;n++)
    printf("%02X", hash[n]);
  printf("\n");
}

void dofile(char *filename, struct stat *statbuf){
  if(alloc_files==num_files){
    alloc_files+=alloc_chunk;
    bloomhashes=realloc(bloomhashes,alloc_files*bloom_bytes);
    filenames=realloc(filenames, alloc_files*sizeof(char *));
  }
  if(!bloomhashes || !filenames){
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  uint8_t *bloomhash=bloomhashes+(num_files*bloom_bytes);
  memset(bloomhash, 0, bloom_bytes);
  int err=hash(filename, bloomhash);
  if(err)
    return;
//  printf("%lu,%s\n", num_files, filename);
//  diag_show_hash(bloomhash);
  filenames[num_files] = malloc(strlen(filename)+1);
  strcpy(filenames[num_files], filename);
  num_files++;
}

int hash(char *filename, uint8_t *bloomhash){
  int fd = open (filename, O_RDONLY);
  if(!fd){
    fprintf(stderr, "Error opening file %s\n", filename);
    return(1);
  }
  struct stat file_info;
  if(fstat (fd, &file_info)){
    fprintf(stderr, "Error statting file %s\n", filename);
    close(fd);
    return(1);
  }  

  if(file_info.st_size==8){
    close(fd);
    return(2);
  }

 uint8_t *mapped_file = mmap (NULL, file_info.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if(!mapped_file){
    fprintf(stderr, "Error mapping file %s\n", filename);
    close(fd);
    return(1);
  }

  uint64_t pos=0;  
  uint64_t filesize=file_info.st_size;
  uint64_t bytes_left=filesize;
  uint32_t bitsset=0;
  do{
    uint32_t chunk_len=get_chunk_len(mapped_file+pos, bytes_left, chunking_threshhold);
    #ifdef diag_mode
      printf("Chunk %08lX %08X\n", pos, chunk_len);
    #endif

    uint64_t fnvhash=make_fnv_from_block(mapped_file+pos, chunk_len);
    uint32_t bucket=fnvhash % (bloom_bytes*8);
    uint32_t bloombyte=bucket>>3;
//    printf("%09lu:%09lu,%06u:%08lX(%u/%u)\n", pos, bytes_left, chunk_len, fnvhash, bucket,bloombyte);
    uint8_t bit=1;
    bit = bit << (bucket&0x07);
//    printf("Change: %u+%02X\n", bloombyte, bit);
    bloomhash[bloombyte] = bloomhash[bloombyte] | bit;
    bitsset++;
    if(bitsset>(bloom_bytes*4)){
      fprintf(stderr, "%s\n  This file is too large to process at current settings. Increase B.\n", filename);
      munmap(mapped_file, file_info.st_size);
      close(fd);
      return(3);
    }
    pos+=chunk_len;
    bytes_left=filesize-pos;
  }while(bytes_left);
  munmap(mapped_file, file_info.st_size);
  close(fd);
  if(bitsset<10){
      fprintf(stderr, "%s\n  This file is too small to process at current settings: Increase T. You may have to increase B in order to handle large files.\n", filename);
      return(4);
  }
  return(0);
}

uint64_t get_file_size(char *filename){
  struct stat statbuf;
  if (stat(filename, &statbuf) == -1){
    fprintf(stderr, "Error statting %s\n", filename);
    return(0);
  }
  return(statbuf.st_size);

}

uint32_t get_chunk_len(uint8_t *data, uint32_t max, uint8_t threshhold){
  uint64_t limit=1;
  limit = limit << threshhold;
  uint32_t n;
  for(n=8;n<max;n++){
    uint64_t hash=make_fnv_from_block(data+n-8, 8);
    if(hash<limit)
      return(n);
    
  }
  return(max);
}

uint64_t make_fnv_from_block(uint8_t *block, uint32_t blocksize){
    uint64_t hash = 14695981039346656037ull;
    for(int n = 0; n < blocksize; n++)
    {
        hash = hash ^ (block[n]);
        hash = hash * 1099511628211; // Magic prime.
    }
    return(hash);
}

void dofolder(char *foldername){
  struct stat statbuf;

  if (stat(foldername, &statbuf) == -1){
    fprintf(stderr, "Error statting %s\n", foldername);
    return;
  }
  if(S_ISREG(statbuf.st_mode)){
    dofile(foldername, &statbuf);
    return;
  }
  if(! S_ISDIR(statbuf.st_mode))
    return; //What is this non-file, non-folder?

  DIR *dp;
  struct dirent *ep;

  dp = opendir (foldername);
  if (dp == NULL){
    fprintf(stderr, "Could not read folder %s\n", foldername);
    return;
  }
  while ((ep = readdir(dp)))
  if(ep->d_name[0]!='.'){
    char *filename=malloc((strlen(foldername)+strlen(ep->d_name)+2)*sizeof(char)); //+2: One for the /, one for the null term.
    strcpy(filename, foldername);
    strcat(filename, "/");
    strcat(filename, ep->d_name);
    dofolder(filename);
    free(filename);
  }
  closedir (dp);
}
