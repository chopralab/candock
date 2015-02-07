#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <json/json.h>

//~ g++ -std=c++0x json.cpp ../jsoncpp-src-0.6.0-rc2/src/lib_json/json_reader.cpp  ../jsoncpp-src-0.6.0-rc2/src/lib_json/json_value.cpp ../jsoncpp-src-0.6.0-rc2/src/lib_json/json_writer.cpp      -o json -I ../jsoncpp-src-0.6.0-rc2/include
//~ Biological Assembiles: ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/
//~ rsync -av --port=33444 ftp.wwpdb.org::ftp_data/biounit .
using namespace std;
int main() {
	//~ Json::Value event;   
	//~ Json::Value vec(Json::arrayValue);
	//~ vec.append(Json::Value(1));
	//~ vec.append(Json::Value(2));
	//~ vec.append(Json::Value(3));
	//~ 
	//~ event["competitors"]["home"]["name"] = "Liverpool";
	//~ event["competitors"]["away"]["code"] = 89223;
	//~ event["competitors"]["away"]["name"] = "Aston Villa";
	//~ event["competitors"]["away"]["code"]=vec;
	//~ event["competitors"]["janez"]["address"]="brilejeva";
	//~ event["competitors"]["janez"]["address"]="cufarjeva";
//~ 
	//~ Json::StyledWriter writer;
	//~ // Make a new JSON document for the configuration. Preserve original comments.
	//~ std::string outputConfig = writer.write( event );
//~ 
	//~ cout << outputConfig << endl;
	//~ std::cout << event << std::endl;

	Json::Value root;   // will contains the root value after parsing.
	Json::Reader reader;
	ifstream infile("../alignments.json");
	string config_doc((istreambuf_iterator<char>(infile)), istreambuf_iterator<char>());
	//~ cout << config_doc << endl;
	bool parsingSuccessful = reader.parse( config_doc, root );
	if ( !parsingSuccessful )
	{
	    // report to the user the failure and their locations in the document.
	    std::cout  << "Failed to parse configuration\n"
	               << reader.getFormattedErrorMessages();
	    return 0;
	}
	for ( int index = 0; index < root.size(); ++index )  { // Iterates over the sequence elements.
		if (root[index]["pdb_id"] == "3e7g") {
			Json::Value vec(Json::arrayValue);
			vec.append("ATP");
			vec.append("NAD");
			vec.append("GDP");
			root[index]["ligands"] = vec;
		}
	}

	Json::FastWriter writer;
	// Make a new JSON document for the configuration. Preserve original comments.
	std::string outputConfig = writer.write( root );
	cout << outputConfig << endl;

	//~ const Json::Value plugins = root;
	//~ for ( int index = 0; index < plugins.size(); ++index )  { // Iterates over the sequence elements.
	   //~ cout << plugins[index]["pdb_id"].asString() << endl;
		//~ const Json::Value alignment = plugins[index]["alignment"];
		//~ for ( int j = 0; j < alignment.size(); ++j ) {
			//~ cout << "\t" << alignment[j]["scores"]["z_score"].asDouble();
			//~ cout << "\t[" << alignment[j]["rotation_matrix"][0][0].asDouble()
				//~ << "," << alignment[j]["rotation_matrix"][0][1].asDouble() 
				//~ << "," << alignment[j]["rotation_matrix"][0][2].asDouble() << "]"
			 //~ << endl;
		//~ }
	//~ }
	//~ cout << root;
	//~ exit(1);
	//~ // Get the value of the member of root named 'encoding', return 'UTF-8' if there is no
	//~ // such member.
	//~ std::string encoding = root.get("encoding", "UTF-8" ).asString();
	//~ // Get the value of the member of root named 'encoding', return a 'null' value if
	//~ // there is no such member.
	//~ const Json::Value plugins = root["plug-ins"];
	//~ for ( int index = 0; index < plugins.size(); ++index )  // Iterates over the sequence elements.
	   //~ cout << plugins[index].asString() << endl;
	//~ setIndentLength( root["indent"].get("length", 3).asInt() );
	//~ setIndentUseSpace( root["indent"].get("use_space", true).asBool() );
	
	// ...
	// At application shutdown to make the new configuration document:
	// Since Json::Value has implicit constructor for all value types, it is not
	// necessary to explicitly construct the Json::Value object:
	//~ root["encoding"] = "janez";
	//~ root["indent"]["length"] = "aska";
	//~ root["indent"]["use_space"] = "nima";
	//~ root["encoding"] = getCurrentEncoding();
	//~ root["indent"]["length"] = getCurrentIndentLength();
	//~ root["indent"]["use_space"] = getCurrentIndentUseSpace();
	
	//~ Json::StyledWriter writer;
	//~ // Make a new JSON document for the configuration. Preserve original comments.
	//~ std::string outputConfig = writer.write( root );
	//~ 
	//~ // You can also use streams.  This will put the contents of any JSON
	//~ // stream at a particular sub-value, if you'd like.
	//~ std::cin >> root["subtree"];
	//~ 
	//~ // And you can write to a stream, using the StyledWriter automatically.
	//~ std::cout << root;
	return 0;
}
