#ifndef JSONREADER_H
#define JSONREADER_H
#include "candock/external/json/json.h"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <math.h>
using namespace std;

class JsonReader {
	Json::Value __root;   // will contains the root value after parsing.
public:
	//~ class iterator : public Json::ValueIterator {};
	typedef Json::ValueIterator iterator;
	iterator begin() { return __root.begin(); }
	iterator end() { return __root.end();}
	void print_JSON_value(const Json::Value&) const;
	void parse_JSON(const string);
	bool print_JSON_tree(const Json::Value &root, const unsigned short depth=0);
	const Json::Value& root() const { return __root; }
	iterator find(const vector<pair<const string, const string>>&); // return node that has all key:value pairs, if not found, return __root.end()
	string output_json() {
		Json::FastWriter writer;
		return writer.write(__root);
	}
};
#endif
