#include "ConsoleProgressBar.hpp"

ConsoleProgressBar::ConsoleProgressBar(std::ostream& stream, size_t repeats, size_t size, const std::string& prefix)
    : _repeats(repeats), _size(size), _current(0), _prefix(prefix), _stream(stream)
{}


void ConsoleProgressBar::display() const
{
    char bar_repr[4096] ;

    double perc = 0. ;
    if(this->_current == this->_repeats)
    {   perc = 100. ; }
    else
    {   perc = static_cast<double>(this->_current) / static_cast<double>(this->_repeats) * 100 ; }
    size_t x = this->_size * this->_current / this->_repeats ;
    std::string bar ;
    for(size_t i=0; i<x; i++) { bar += "=" ; }
    for(size_t i=x; i<this->_size; i++) { bar += "." ; }

    sprintf(bar_repr, "%s : progress [%s] %.2f %%\r", this->_prefix.c_str(), bar.c_str(), perc) ;
    this->_stream << bar_repr ;
    this->_stream.flush() ;
}


void ConsoleProgressBar::update()
{   this->display() ;
    if(this->_current < this->_repeats)
    {   this->_current++ ; }
}
