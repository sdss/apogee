#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 

void error(char *msg)
{
    perror(msg);
    exit(0);
}

int main(int argc, char *argv[])
{
    int sockfd, portno=1050, base=2,n,i,retval,ntot;

    struct sockaddr_in serv_addr;
    struct hostent *server;

#define MAXLEN 1000000
    char buffer[MAXLEN];
    if (argc < 3) {
       fprintf(stderr,"usage %s hostname portno\n", argv[0]);
       exit(0);
    }
    server=gethostbyname(argv[1]);
    sscanf(argv[2],"%d",&portno);
    base=3;
fprintf(stderr,"port: %d %d %d\n",portno,argc,base);

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) 
        error("ERROR opening socket");
    if (server == NULL) {
        fprintf(stderr,"ERROR, no such host\n");
        exit(0);
    }
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr, 
         (char *)&serv_addr.sin_addr.s_addr,
         server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(sockfd,(struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0) 
        error("ERROR connecting");
    bzero(buffer,MAXLEN);
    //printf("Please enter the message: ");
    //fgets(buffer,255,stdin);
    for (i=base; i<argc ; i++) sprintf(buffer+strlen(buffer),"%s ",argv[i]);
    n = write(sockfd,buffer,strlen(buffer));
    if (n < 0) 
         error("ERROR writing to socket");
    bzero(buffer,MAXLEN);
    n=-1;
    ntot=0;
    while (n != 0) {
     n = read(sockfd,buffer+ntot,MAXLEN-1);
     ntot+=n;
    }
    fprintf(stderr,"%s\n",buffer);
    sscanf(buffer,"%d",&retval);

    if (n < 0) 
         error("ERROR reading from socket");
    close(sockfd);
    return(retval);
}
