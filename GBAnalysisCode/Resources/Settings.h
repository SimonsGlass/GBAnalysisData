/*
 *  Settings.h
 *  Percolation
 *
 *  Created by Samuel Schoenholz on 6/4/10.
 *  Copyright 2010 Swarthmore College. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <map>

#include "Value.h"

using namespace std;
//#define DEBUGTAS

//Class to read settings from a file
class Settings
{
protected:
	
	map<string,Value> ItemMap;

public:
	Settings() {}
	
	//add a settings item
	void AddItem(string,Value);
	//find a settings item
	Value operator[](string token);

    //reads settings from command line arguments
    void Read(int argc, char *argv[]);
    
	//Read settings from a stream.
	friend istream &operator>>(istream &in, Settings &settings)
	{
		map<string,Value>::iterator it;
		string token;
		
		//while either we have not yet reached the end of the file or we have not found a token "end"
		while(!in.eof() && token.compare("ENDFILE")!=0)
		{
			in >> token;
#ifdef DEBUGTAS
                        cout << "We read in a new token: " << token << endl; 
#endif

                        // Skip words between the QQ tokens
                        if (token.compare("QQ") == 0) {
                                do in >> token; 
                                while (token.compare("QQ")!=0 && !in.eof() && token.compare("ENDFILE")!=0);
                        } else {
#ifdef DEBUGTAS
                            cout << "Settings file token " << token << endl;
#endif
    			    //If the current token maps to an Item then read the Item into the map
			    it = settings.ItemMap.find(token);
			    if(it!=settings.ItemMap.end()){
				Value itm = (*it).second;
				in >> itm;
                                cout << "Settings file token " << token;
				cout << " matched as an Item and is given value " << itm;
				settings.ItemMap[token] = itm;
			        cout << endl;
			    }
                        }
		}
		return in;
	}
	
	//Write settings to a stream
	friend ostream &operator<<(ostream &out, Settings &settings)
	{
		map<string,Value>::iterator it;
		for(it = settings.ItemMap.begin();it!=settings.ItemMap.end();it++)
		{
			out << (*it).first << " " << (*it).second << endl;
		}
		out << "end\n";
		return out;
		
	}
};

//reads settings from command line arguments
void Settings::Read(int argc, char *argv[])
{
    //for(int i = 0;i<argc;i++)
      //  argv[i];
}

//adds an item to the map
void Settings::AddItem(string token,Value item)
{
	ItemMap[token] = item;
}

//finds an item and returns it
Value Settings::operator[](string token)
{
	//Find an item using a token and return it, otherwise return null.
	map<string,Value>::iterator it;
	it = ItemMap.find(token);
	if(it!=ItemMap.end())
		return (*it).second;
        cout << "We did not find an Item " << token << endl;
	return Value();
}


