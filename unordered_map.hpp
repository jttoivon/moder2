/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

    MODER is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MODER is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/
#include <boost/unordered_map.hpp>


// This version doesn't add default value to map
// if key is not found in map
template <typename K, typename V, typename H = boost::hash<K> >
class my_unordered_map
  : public boost::unordered_map<K, V, H>
{
public:

  my_unordered_map() : boost::unordered_map<K, V, H>() {}

  V
  operator[](K k) const
  {
    typename boost::unordered_map<K, V, H>::const_iterator i = this->find(k);
    if (i == boost::unordered_map<K, V, H>::end())
      return V();
    else
      return i->second;
  }

  V&
  operator[](K k)
  {
    return boost::unordered_map<K, V, H>::operator[](k);
  }

};


