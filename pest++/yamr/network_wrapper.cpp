#include "network_wrapper.h"
#include "network_package.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>

#ifdef OS_LINUX
  #include <arpa/inet.h>
  #include <unistd.h>
#endif

using namespace std;

void w_init()
{
#ifdef OS_WIN
	WSADATA wsaData;

	if (WSAStartup(MAKEWORD(2,0),  &wsaData)!=0)
	{
		cerr << "WSAStartup failed " << w_get_error_msg() << endl;
	}
#endif
}

int w_close(int sockfd) 
{
  int n;
    #ifdef OS_WIN
	shutdown(sockfd, SD_BOTH);
	if ((n = closesocket(sockfd)) != 0)
	{
		cerr << "error closing socket: "  << w_get_error_msg() << endl;  
	}
	return n;
    #endif
    #ifdef OS_LINUX
	n = shutdown(sockfd, 2);
    #endif
    return n;
}
void w_cleanup()
{
   #ifdef OS_WIN
	WSACleanup();
   #endif
}


int w_getaddrinfo(const char *node, const char *service,
			  const struct addrinfo *hints, struct addrinfo **res)
{
	int status;
	if ((status = getaddrinfo(node, service, hints, res)) !=0)
	{
		cerr << "getaddrinfo error: " << gai_strerror(status) << endl;
	}
	return status;
}

vector<string> w_getnameinfo_vec(int sockfd, int flags)
{
	int err;
	vector<string> name_info;
	char host[INET6_ADDRSTRLEN];
	char port[INET6_ADDRSTRLEN];
	struct sockaddr_storage addr;
	socklen_t addr_len = sizeof addr;
	err = getpeername(sockfd, (struct sockaddr*) &addr, &addr_len);
	err = getnameinfo((struct sockaddr*) &addr, addr_len, host, sizeof host, port, sizeof port, flags); 
	name_info.push_back(host);
	name_info.push_back(port);
	return name_info;
}

int w_socket(int domain, int type, int protocol)
{
	int sockfd = socket(domain, type, protocol);
	if (sockfd < 0) {
		cerr << "socket error: "  << w_get_error_msg() << endl;
	}
	return sockfd;
}

int w_connect(int sockfd, struct sockaddr *serv_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=connect(sockfd, serv_addr, addrlen)) == -1 )
	{
		cerr << "connect error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_bind(int sockfd, struct sockaddr *my_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=bind(sockfd, my_addr, addrlen)) == -1 )
	{
		cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_accept(int sockfd, struct sockaddr *addr, socklen_t *addrlen)
{
	int n=0;
	if ((n=accept(sockfd, addr, addrlen)) == -1) 
	{
		cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}


int w_listen(int sockfd, int backlog)
{
	int n;
	if ((n = listen(sockfd, backlog)) == -1)
	{
		cerr << "listen error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_recv(int sockfd, char *buf, size_t len, int flags)
{
	int n;
	n = recv(sockfd, buf, len, flags);
	if (n < 0){
		cerr << "recv error: "  << w_get_error_msg() << endl;
	}
	return n;
}
int w_send(int sockfd, char *buf, size_t len, int flags)
{
	int n;
	n = send(sockfd, buf, len, flags);
	if (n < 0){
		cerr << "send error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_sendall(int sockfd, char *buf, unsigned long *len)
{
	unsigned long total = 0; // how many bytes we've sent
	unsigned long bytesleft = *len; // how many we have left to send
	int n;
	while(total < *len) {
		n = send(sockfd, buf+total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually sent here
		if (n < 0){
		cerr << "w_sendall error: " << n << endl;
	}
	if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}


int w_recvall(int sockfd, char *buf, unsigned long *len)
{
	unsigned long total = 0; // how many bytes we've received
	unsigned long bytesleft = *len; // how many we have left to receive
	int n;
	while(total < *len) {
		n = recv(sockfd, buf+total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually received here
	if (n < 0){
		cerr << "w_recvall error: " << n << endl;
	}
	if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}



//addrinfo* w_bind_first_avl(addrinfo *servinfo)
//{
//	// !!!!!!!!!!!!!!!THIS NEEDS REVIEW AND TESTING !!!!!!!!!!!!!!!!!!!!!
//	//
//	// loop through all the results and bind to the first we can
//	struct addrinfo *p;
//	int sockfd;
//	char yes = '1';
//	for(p = servinfo; p != NULL; p = p->ai_next) 
//	{
//		if ((sockfd = socket(p->ai_family, p->ai_socktype,
//		p->ai_protocol)) == -1) 
//		{
//			perror("server: socket");
//			continue;
//		}
//		if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes,
//		sizeof(int)) == -1)
//		{
//			perror("setsockopt");
//			return NULL;
//		}
//		if (bind(sockfd, p->ai_addr, p->ai_addrlen) == -1) 
//		{
//			w_close(sockfd);
//			perror("server: bind");
//			continue;
//		}
//	break;
//	}
//	if (p == NULL) {
//		cerr << "server: failed to bind" << endl;
//	}
//	return p;
//}

void w_print_servinfo(addrinfo *res, ostream &fout)
{
	struct addrinfo *p;
	fout << "IP addresses:" << endl;
	for(p = res;p != NULL; p = p->ai_next) 
	{
		void *addr;
		char ipstr[INET6_ADDRSTRLEN];
		string ipver;
		// get the pointer to the address itself,
		// different fields in IPv4 and IPv6:
		if (p->ai_family == AF_INET) { // IPv4
			struct sockaddr_in *ipv4 = (struct sockaddr_in *)p->ai_addr;
			addr = &(ipv4->sin_addr);
			ipver = "IPv4";
		}
		else { // IPv6
			struct sockaddr_in6 *ipv6 = (struct sockaddr_in6 *)p->ai_addr;
			addr = &(ipv6->sin6_addr);
			ipver = "IPv6";
		}
	// convert the IP to a string and print it:
	inet_ntop(p->ai_family, addr, ipstr, sizeof ipstr);
	fout << " " << ipver <<": " << ipstr << endl;
	}
}


int w_select(int numfds, fd_set *readfds, fd_set *writefds,
		   fd_set *exceptfds, struct timeval *timeout)
{
	int n;
	if ((n=select(numfds, readfds, writefds, exceptfds, timeout)) == -1)
	{
		cerr << "select error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_memcpy_s(void *dest, size_t numberOfElements, const void *src, size_t count)
{
	int err;
        #ifdef OS_WIN
	err = memcpy_s(dest, numberOfElements, src, count);
        #endif
        #ifdef OS_LINUX
	memcpy(dest, src, count);
        #endif
	if (err) {
		cerr << "Error executing memcpy" << endl;
	}
	return err;
}

string w_get_error_msg()
{
	stringstream err_msg;
	#ifdef OS_WIN
	int err_no = WSAGetLastError();
	LPVOID lpMsgBuf;
	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER |
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		err_no,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
		(LPTSTR) &lpMsgBuf,
		0,
		NULL);
	// Process any inserts in lpMsgBuf.
	// ..
	err_msg << "(tcp/ip error id = " << err_no << ") " << (LPCTSTR)lpMsgBuf ;
	// Free the buffer.
	LocalFree( lpMsgBuf );
	#endif
	#ifdef OS_LINUX
	err_msg << "(tcp/ip error id = " << errno << ")";
	#endif
	return err_msg.str();
}

void w_sleep(int millisec)
{
   #ifdef OS_WIN
	Sleep(millisec);
   #endif
   #ifdef OS_LINUX
	sleep(millisec / 1000);
   #endif
}
