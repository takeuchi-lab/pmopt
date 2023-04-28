#include "pmopt/databases/sequence_database.hpp"

#include <cassert>
#include <iostream>
#include <sstream>

using namespace PMOpt;

void SequenceDatabase::push(std::string && buffer)
{
    if (_format == FORMAT_STRING)
        push_string(std::move(buffer));
    else if (_format == FORMAT_SEQUENCE)
        push_sequence(std::move(buffer));
}


void SequenceDatabase::push_string(std::string && buffer)
{
    Sequence sequence;

    for (auto c : buffer)
    {
        auto item = std::string{c};
        if (_item2id.find(item) == _item2id.end())
        {
            // found new item
            _item2id[item] = _items.size();
            _items.push_back(item);
        }
        sequence.push_back(_item2id[item]);
    }

    emplace_back(std::move(sequence));
}


void SequenceDatabase::push_sequence(std::string && buffer)
{
    Sequence sequence;

    std::stringstream ss(std::move(buffer));
    while (std::getline(ss, buffer, ' '))
    {
        if (_item2id.find(buffer) == _item2id.end())
        {
            // found new item
            _item2id[buffer] = _items.size();
            _items.push_back(buffer);
        }
        sequence.push_back(_item2id[buffer]);
    }

    emplace_back(std::move(sequence));
}

std::string SequenceDatabase::struct_to_str(const Sequence & structure) const noexcept
{
    if (_format == FORMAT_STRING)
        return struct_to_str_string(structure);
    else
        return struct_to_str_sequence(structure);
}


std::string SequenceDatabase::struct_to_str_sequence(const Sequence & structure) const noexcept
{
    std::stringstream ss;

    bool is_first = true;
    for (auto item : structure)
    {
        assert(item != ROOT);

        if (!is_first)
            ss << " ";
        ss << _items[item];
        is_first = false;
    }

    return ss.str();
}

std::string SequenceDatabase::struct_to_str_string(const Sequence & structure) const noexcept
{
    std::stringstream ss;

    for (auto item : structure)
    {
        assert(item != ROOT);
        ss << _items[item];
    }

    return ss.str();
}