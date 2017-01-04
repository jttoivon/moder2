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


