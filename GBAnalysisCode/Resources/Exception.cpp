#include "Exception.h"

namespace LiuJamming
{

CException::CException()
{
	Function = "";
	Error = "";
}

CException::CException(string _Function, string _Error)
{
	Function = _Function;
	Error = _Error;
}
	
const string &CException::GetFunction() const
{
	return Function;
}

const string &CException::GetError() const
{
	return Error;
}
	
std::ostream &operator<<(std::ostream &out, const CException &exception)
{
	out << "In function: " << exception.GetFunction() << " encountered the error message:\n\"" << exception.GetError() << "\"\n";
	return out;
}

}