#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <sstream>
#include <algorithm>
#include "jsonreader.hpp"
#include "helper/error.hpp"
#include "helper/inout.hpp"
#include <lib_json/json_reader.cpp>
#include <lib_json/json_value.cpp>
#include <lib_json/json_writer.cpp>

JsonReader::iterator JsonReader::find(const vector<pair<const string, const string> > &kv) {
	for(Json::ValueIterator itr = __root.begin() ; itr != __root.end() ; itr++ ) {
		unsigned int i = 0;
		for (auto &k : kv) {
			const string key = k.first;
			const string value = k.second;
			if (value == (*itr)[key].asString()) i++;
		}
		if (i == kv.size()) return itr;
	}
	return this->end();
}

void JsonReader::print_JSON_value(const Json::Value &val) const {
    if( val.isString() ) {
        cout << val.asString(); 
    } else if( val.isBool() ) {
        cout << val.asBool(); 
    } else if( val.isInt() ) {
        cout << val.asInt(); 
    } else if( val.isUInt() ) {
        cout << val.asUInt(); 
    } else if( val.isDouble() ) {
        cout << val.asDouble(); 
    }
    else {
        cout << "unknown type=[" << val.type() << "]"; 
    }
}

bool JsonReader::print_JSON_tree(const Json::Value &root, const unsigned short depth) {
    cout << " {type=["<< root.type() << "], size=" << root.size() << "}"; 
    if( root.size() > 0 ) {
        cout<<endl;
        for(Json::ValueConstIterator itr = root.begin() ; itr != root.end() ; itr++ ) {
            // Print depth. 
            for( int tab = 0 ; tab < depth; tab++) {
               cout << "-";
            }
            cout << " subvalue(";
            print_JSON_value(itr.key());
            cout << ") -";
            print_JSON_tree(*itr, depth + 1); 
        }
        return true;
    } else {
        cout << " ";
        print_JSON_value(root);
        cout << "\n";
    }
    return true;
}

void JsonReader::parse_JSON(const string JSON_file) {
	if (boost::regex_search(JSON_file, boost::regex(".json$"))) {
		Json::Reader reader;
		vector<string> vec;
		inout::Inout::read_file(JSON_file, vec);
		ostringstream ss;
		copy(vec.begin(), vec.end(), ostream_iterator<string>(ss, ""));
		if (!reader.parse(ss.str(), __root)) {
			throw Error("Failed to parse JSON file " + reader.getFormattedErrorMessages());
		}
	} 
	else {
		throw Error("die : need a .json file");
	}	
}

