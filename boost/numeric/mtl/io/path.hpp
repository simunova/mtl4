// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef MTL_IO_PATH_INCLUDE
#define MTL_IO_PATH_INCLUDE

namespace mtl { namespace io {

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
    const static char delim = '\\';
#else
    const static char delim = '/';
#endif

/// Join the directory and file
/** It concatenates both with the os-specific slash unless directory is empty, then only file is returned **/
std::string inline join(std::string directory, std::string file)
{
    return directory.empty() ? file : directory + delim + file;
}

/// Directory name in s that is everything before last slash
std::string inline directory_name(std::string s)
{
    for (int i= int(s.size()) - 1; i >= 0; i--)
	if (s[i] == delim)
	    return s.substr(0, i);
    return std::string();
}

/// File name in s that is everything after last slash
std::string inline file_name(std::string s)
{
    for (int i= int(s.size()) - 1; i >= 0; i--)
	if (s[i] == delim)
	    return s.substr(i + 1);
    return s;
}




}} // namespace mtl::io

#endif // MTL_IO_PATH_INCLUDE
