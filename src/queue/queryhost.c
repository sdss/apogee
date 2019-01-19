/* host program to keep track of job IDs that are running and/or done 
     input should be of  form:   command id set
      where command is stop, status, done, running
      set is 0 to query status, 1 to set status
      id is character job ID up to 64 chars
     returns 0 on succesful query/set, else 1
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <errno.h>
#define MAXPROC 35000
int main(int argc, char *argv[])
{
  char done[MAXPROC][64], running[MAXPROC][64], id[64];
  int ndone=0, nrunning=0, set,retval,i;
  char query[64],command[80], retstring[1000000];
  //FILE *socket;

  int sockfd, newsockfd, portno, clilen;
  char buffer[256];
  struct sockaddr_in serv_addr, cli_addr;
  int n, port=1050, irun, idone, verbose=0;

  if (argc < 1) {
      fprintf(stderr,"usage %s portno\n", argv[0]);
      exit(0);
   }
   if (argc==2) sscanf(argv[1],"%d",&port); 

  // set up socket host
  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) perror("ERROR opening socket");
  bzero((char *) &serv_addr, sizeof(serv_addr));
  portno = port;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(portno);
  if (bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
              error("ERROR on binding");
  listen(sockfd,5);

  // continue to wait for status requests / notifications
  while (1) {
    clilen = sizeof(cli_addr);
    newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr, &clilen);
    if (newsockfd < 0) perror("ERROR on accept");
    bzero(retstring,1000000);
    bzero(buffer,256);
    bzero(command,80);
    n = read(newsockfd,command,80);
    //socket=fopen("query","r");
    //fgets(command,80,socket);
    sscanf(command,"%s%s%d",query,id,&set);
    if (verbose) {
      fprintf(stderr,"received %s %s %d\n",query,id,set);
      fprintf(stderr,"%d %d\n",nrunning,ndone);
    }

    if (strcmp(query,"stop")==0) exit(0);
    
    // status command simply prints what is done and what is running
    if (strcmp(query,"status")==0) {
      if (verbose) fprintf(stderr,"%s\n",query);
      sprintf(retstring,"%d",0);
      sprintf(retstring+strlen(retstring),"\n done: %d",ndone);
      for (i=0; i<ndone; i++) sprintf(retstring+strlen(retstring)," %s",done[i]);
      sprintf(retstring+strlen(retstring),"\n running: %d",nrunning);
      for (i=0; i<nrunning; i++) if (strcmp(running[i],"0") !=0 ) sprintf(retstring+strlen(retstring)," %s",running[i]);
      if (verbose) fprintf(stderr,"strlen: %d\n",strlen(retstring));
      write(newsockfd,retstring,strlen(retstring));
      retval=-1;
    }

    // done command
    if (strcmp(query,"done") == 0) {
      if (verbose) fprintf(stderr,"%s %d\n",query,set);
      if (set==1) {
        irun=-1;
        for (i=0; i<nrunning; i++) if (strcmp(running[i],id) == 0) irun=i;
        if (irun >= 0) {
          for (i=irun; i<nrunning-1; i++) strcpy(running[i],running[i+1]);
          nrunning--;
          if (verbose) fprintf(stderr,"done set: %d\n",nrunning);
        }
        strcpy(done[ndone++],id);
        retval=0;
      } else {
        retval=1;
        for (i=0; i<ndone; i++) if (strcmp(done[i],id) == 0) retval=0;
      }
    }

    // running command
    if (strcmp(query,"running") == 0) {
      if (verbose) fprintf(stderr,"%s %d\n",query, set);
      if (set==1) {
        // if we get a request to run, confirm that this ID is not already done or running
        retval=0;
        for (i=0; i<ndone; i++) if (strcmp(done[i],id) == 0) retval=1;
        for (i=0; i<nrunning; i++) if (strcmp(running[i],id) == 0) retval=1;
        if (retval==0) strcpy(running[nrunning++],id);
      } else {
        retval=1;
        for (i=0; i<nrunning; i++) {
         if (strcmp(running[i],id) == 0) retval=0;
        }
      }
    }
    if (strcmp(query,"clear") == 0) {
      if (verbose) fprintf(stderr,"%s %d\n",query,set);
      if (strcmp(id,"all") ==0) {
        nrunning = 0;
        ndone = 0;
      } else {
        irun=-1;
        for (i=0; i<nrunning; i++) if (strcmp(running[i],id) == 0) irun=i;
        if (irun >= 0) {
          for (i=irun; i<nrunning-1; i++) strcpy(running[i],running[i+1]);
          nrunning--;
          if (verbose) fprintf(stderr,"clear running : %d\n",nrunning);
        }
        idone=-1;
        for (i=0; i<ndone; i++) if (strcmp(done[i],id) == 0) idone=i;
        if (idone >= 0) {
          for (i=idone; i<ndone-1; i++) strcpy(done[i],done[i+1]);
          ndone--;
          if (verbose) fprintf(stderr,"clear done : %d\n",ndone);
        }
      }
      retval=0;
    }
   
    // send the return value back to the client
    if (retval>=0) {
      sprintf(retstring,"%d",retval);
      write(newsockfd,retstring,strlen(retstring));
    }
    close(newsockfd);

    //fclose(socket);
  }
}
