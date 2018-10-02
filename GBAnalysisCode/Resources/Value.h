/*
 *  Value.h
 *  Percolation
 *
 *  Created by Samuel Schoenholz on 6/7/10.
 *  Copyright 2010 Swarthmore College. All rights reserved.
 *
 */

#ifndef VALUE
#define VALUE

#include <map>
#include <string>
#include <typeinfo>

/*
This code introduces typeless, invisible variables.
The code is constructed from the article
Chameleon Objects, or how to write a generic, type safe wrapper class
Volker Simonis1 
Wilhelm-Schickard-Institut für Informatik 
Universität Tübingen, 72076 Tübingen, Germany 
E-mail: simonis@informatik.uni-tuebingen.de
 
29 November 1999
*/

//class that handles casting between incompatible types
class Incompatible_Type_Exception
{
private:
	string errorType;
public:
	Incompatible_Type_Exception(const string &s) {errorType = s;}
	string getError() const {return errorType;}
};

template <class T>
istream &operator>>(istream &in,T *v)
{
	return in >> *v;
}

//class that seamlessly handles arbitrary types
class Value
{
protected:
	enum Action {GET,SET,DELETE};
	//primary method that adds items to a static lookup table
	template<class T> T& value(T t = T(),Action = GET) throw (Incompatible_Type_Exception&);
	typedef void(Value::*FuncPointer)(void);
	typedef void(Value::*ClonePointer)(Value*) const;
	typedef ostream&(Value::*PrintPointer)(ostream&) const;
	typedef istream&(Value::*ReadPointer)(istream&) const;
	typedef Value(Value::*DereferencePointer)(void) const;
	ClonePointer fp_CLONE;								//pointer to the cloning function
	FuncPointer fp_DELETE;								//pointer to the deletion function
	PrintPointer fp_PRINT;								//pointer to the printer function
	ReadPointer fp_READ;								//pointer to the writing function
	DereferencePointer fp_DEREFERENCE;
	
public:	
	enum Type {POINTER};

	Value() : fp_DELETE(NULL), fp_CLONE(NULL) , fp_PRINT(NULL), fp_READ(NULL), fp_DEREFERENCE(NULL) {}		//default constructor
	template<class T> Value(const T& t) : fp_DELETE(NULL), fp_DEREFERENCE(NULL) {value(t,SET);}	//generic constructor
	template<class T> Value(T *t,Type) : fp_DELETE(NULL) {value<T*>(t,SET); fp_DEREFERENCE = &Value::DereferenceValue<T>;}
	
	inline Value(const Value &);						//copy constructor
	~Value() {if(fp_DELETE!=NULL) (this->*fp_DELETE)();}//generic destructor
	template<class T> operator T() const throw (Incompatible_Type_Exception&) {	//generic casting function
		return const_cast<Value*>(this)->template value<T>();
	}
	template <class T> T& operator=(const T &t) {return value(t,SET);}
	Value &operator=(const Value&);						//clone function
	template <class T> inline void DeleteValue();		//helper function to delete value
	template <class T> void CloneValue(Value *) const;	//generic function to clone an object
	
	template <class T> void Pointer(T *t) {value<T*>(t,SET);fp_DEREFERENCE = &Value::DereferenceValue<T>;}
	Value operator*();
	template <class T> Value DereferenceValue() const;

	//printing functions
	template <class T> ostream &PrintValue(ostream &out) const{
		return out << const_cast<Value*>(this)->template value<T>();}
	friend ostream &operator <<(ostream &out, const Value &e){
		if(e.fp_PRINT) return (e.*e.fp_PRINT)(out);
		else return (out << "nil");
	}
	
	//writing functions
	template <class T>istream &ReadValue(istream &in) const{
		T t = const_cast<Value*>(this)->template value<T>(); in >> t; const_cast<Value*>(this)->template value<T>(t,SET);return in;
	}
	friend istream &operator>>(istream &in, const Value &e){
		if(e.fp_READ) return (e.*e.fp_READ)(in);
		else return in;
	}
};

//copy creation operator
inline Value::Value(const Value &val)	: fp_DELETE(NULL), fp_PRINT(NULL), fp_READ(NULL), fp_DEREFERENCE(NULL)
{
	if(val.fp_CLONE!=NULL)
		(val.*val.fp_CLONE)(this);
	else fp_CLONE = NULL;
}

//implementation of value uses a static lookup table for each datatype to provide functionality
template <class T>
T &Value::value(T t,Action action) throw (Incompatible_Type_Exception&)
{
	static std::map<Value*, T> values;
	switch(action)
	{
		case SET: //sets an object in the map for the correct type of object
		{
			values[this] = t;
			fp_READ = &Value::ReadValue<T>;
			fp_PRINT = &Value::PrintValue<T>;
			fp_CLONE = &Value::CloneValue<T>;//sets the clone pointer to have the correct type
			if(fp_DELETE == NULL) fp_DELETE = &Value::DeleteValue<T>;	//set fp_DELETE to point to this types function pointer
			else if (fp_DELETE!=(FuncPointer)&Value::DeleteValue<T>)
			{
				(this->*fp_DELETE)();		//otherwise delete the type
				fp_DELETE = &Value::DeleteValue<T>;//and remember the deletion function for the correct type
			}
			return t;
		}
		case GET: //if the object of this type exists then get it otherwise throw and error
		{
			if(values.count(this))
				return values[this];
			else 
				throw Incompatible_Type_Exception(typeid(T).name());
		}
		case DELETE: //only called by the destructor
		{
			if(values.count(this)) values.erase(this);
			return t;
		}
		default:
			return t;
	}
}

//function to delete a value of type T
template<class T>
inline void Value::DeleteValue()
{
	static T t;
	value(t,DELETE);
}

//function to clone Value.
//fp_CLONE holds a pointer to the function of type val.
Value &Value::operator=(const Value& val)
{
	if(this != &val&&val.fp_CLONE!=NULL){
		(val.*val.fp_CLONE)(this);
		return *this;
	}
}

//function to actually perform the cloning of an object
template <class T>
void Value::CloneValue(Value *val) const
{
	val->template value<T>(const_cast<Value*>(this)->template value<T>(),SET);
	if(fp_DEREFERENCE){val->fp_DEREFERENCE = fp_DEREFERENCE;}
	
}

//operator to give a value that is one time dereferenced from the value
Value Value::operator*()
{
	if(fp_DEREFERENCE)
		 return (this->*fp_DEREFERENCE)();
}

//function to return a dereferenced copy of this value
template<class T>
Value Value::DereferenceValue() const
{
	return Value(*(const_cast<Value*>(this)->template value<T*>()));
}

#define INTEGER 0
#define STRING string()
#define FLOAT 0.0f
#define DOUBLE 0.0

#endif
