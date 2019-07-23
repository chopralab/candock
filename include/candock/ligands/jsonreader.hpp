/* This is jsonreader.hpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef JSONREADER_H
#define JSONREADER_H
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include "candock/external/json/json.h"

namespace candock {

class JsonReader {
    Json::Value __root;  // will contains the root value after parsing.
   public:
    //~ class iterator : public Json::ValueIterator {};
    typedef Json::ValueIterator iterator;
    iterator begin() { return __root.begin(); }
    iterator end() { return __root.end(); }
    void print_JSON_value(const Json::Value&) const;
    void parse_JSON(const std::string);
    bool print_JSON_tree(const Json::Value& root,
                         const unsigned short depth = 0);
    const Json::Value& root() const { return __root; }
    iterator find(
        const std::vector<
            std::pair<const std::string, const std::string>>&);  // return node
                                                                 // that has all
                                                                 // key:value
                                                                 // pairs, if
                                                                 // not found,
                                                                 // return
                                                                 // __root.end()
    std::string output_json() {
        Json::FastWriter writer;
        return writer.write(__root);
    }
};
}

#endif
