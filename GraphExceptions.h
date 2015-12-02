#include <iostream>
#include <exception>

using namespace std;

class GraphException: public exception
{
    virtual const char* what() const throw(){
	return "GraphException thrown";
    }

};
