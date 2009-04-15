/*
 * Copyright (C) 2005 Universitat d'Alacant / Universidad de Alicante
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */
//#include <apertium/hmm2.h>
//#include "hmm2.h"
//#include <apertium/tagger_data_trigram.h>
#include "tagger_data_trigram.h"
#include <lttoolbox/compression.h>
#include <apertium/endian_double_util.h>
#include <apertium/string_utils.h>

using namespace Apertium;

void
TaggerDataTrigram::copy(TaggerDataTrigram const &o)
{
  open_class = o.open_class;
  forbid_rules = o.forbid_rules;
  tag_index = o.tag_index;
  array_tags = o.array_tags;
  enforce_rules = o.enforce_rules;
  prefer_rules = o.prefer_rules;
  constants = o.constants;
  output = o.output;  
  this->setProbabilities(o.N, o.M, o.a, o.b);
  plist = o.plist;
}

void
TaggerDataTrigram::destroy()
{
  if(a != NULL)
  {
    for(int i = 0; i != N; i++)
    {
    	for(int j = 0; j != N; j++)
    	{
      		delete [] a[i][j];
    	}
    }
    delete [] a;
  }
  a = NULL;

  if(b != NULL)
  {
    for(int i = 0; i != N; i++)
    {
    	for(int j = 0; j != N; j++)
	{
      		delete [] b[i][j];
	}
    }
    delete [] b;
  }
  b = NULL;
  N = 0;
  M = 0;
}

TaggerDataTrigram::TaggerDataTrigram()
{
  a = NULL;
  b = NULL;
  N = 0;
  M = 0;
}

TaggerDataTrigram::~TaggerDataTrigram()
{
  destroy();
}

TaggerDataTrigram::TaggerDataTrigram(TaggerDataTrigram const &o)
{
  a = NULL;
  b = NULL;
  N = 0;
  M = 0;
  copy(o);
}

TaggerDataTrigram &
TaggerDataTrigram::operator =(TaggerDataTrigram const &o)
{
  if(this != &o)
  {
    destroy();
    copy(o);
  }
  return *this;
}

set<TTag> &
TaggerDataTrigram::getOpenClass()
{
  return open_class;
}

void
TaggerDataTrigram::setOpenClass(set<TTag> const &oc)
{
  open_class = oc;
}

vector<TForbidRule> &
TaggerDataTrigram::getForbidRules()
{
  return forbid_rules;
}

void
TaggerDataTrigram::setForbidRules(vector<TForbidRule> &fr)
{
  forbid_rules = fr;
}  

map<wstring, TTag, Ltstr> &
TaggerDataTrigram::getTagIndex()
{
  return tag_index;
}

void
TaggerDataTrigram::setTagIndex(map<wstring, TTag, Ltstr> const &ti)
{
  tag_index = ti;
}
  
vector<wstring> &
TaggerDataTrigram::getArrayTags()
{
  return array_tags;
}

void
TaggerDataTrigram::setArrayTags(vector<wstring> const &at)
{
  array_tags = at;
}

vector<TEnforceAfterRule> &
TaggerDataTrigram::getEnforceRules()
{
  return enforce_rules;
}


void
TaggerDataTrigram::setEnforceRules(vector<TEnforceAfterRule> const &tear)
{
  enforce_rules = tear;
}

vector<wstring> &
TaggerDataTrigram::getPreferRules()
{
  return prefer_rules;
}

void
TaggerDataTrigram::setPreferRules(vector<wstring> const &pr)
{
  prefer_rules = pr;
}

vector<wstring> &
TaggerDataTrigram::getDiscardRules()
{
  return discard;
}

void
TaggerDataTrigram::setDiscardRules(vector<wstring> const &v)
{
  discard = v;
}

ConstantManager &
TaggerDataTrigram::getConstants()
{
  return constants;
}

void
TaggerDataTrigram::setConstants(ConstantManager const &c)
{  
  constants = c;
}

Collection &
TaggerDataTrigram::getOutput()
{
  return output;
}

void
TaggerDataTrigram::setOutput(Collection const &c)
{
  output = c;
}
  
void
TaggerDataTrigram::setProbabilities(int const myN, int const myM, 
                             double ***myA, double ***myB)
{
  int i, i1, i2, j;
  this->destroy();
  N = myN;
  M = myM;
  
  if(N != 0 && M != 0)
  {
    // NxNxN matrix
    a = new double ** [N];
    for(i = 0; i != N; i++)
    {
	    a[i] = new double *[N];
    	for(i1 = 0; i1 != N; i1++)
      {
        a[i][i1] = new double[N];
	if(myA != NULL)
	{
	  for(j = 0; j != N; j++)
	  { 
	    a[i][i1][j] = myA[i][i1][j];
	  }
	}
      }
    }
  
    // NxNxM matrix
    b = new double ** [N];
    for(i = 0; i != N; i++)
    {
      b[i] = new double *[N];
      for(i2 = 0; i2 != N; i2++)
      {
        b[i][i2] = new double[M];
        if(myB != NULL)
        {
          for(j = 0; j != M; j++)
          {
            b[i][i2][j] = myB[i][i2][j];
       	  }
        }
      }
    }
  }
  else
  {
    a = NULL;
    b = NULL;
  }  
}

double *** 
TaggerDataTrigram::getA()
{
  return a;
}

double ***
TaggerDataTrigram::getB()
{
  return b;
}

int 
TaggerDataTrigram::getN()
{  
  return N;
}

int
TaggerDataTrigram::getM()
{
  return M;
}

PatternList &
TaggerDataTrigram::getPatternList()
{
  return plist;
}

void
TaggerDataTrigram::setPatternList(PatternList const &pl)
{
  plist = pl;
}

void
TaggerDataTrigram::read(FILE *in)
{
  destroy();

  // open_class
  int val = 0;
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    val += Compression::multibyte_read(in);
    open_class.insert(val);
  }
  wcerr<<L"TaggerDataTrigram::read open_class done\n";
  
  // forbid_rules
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    TForbidRule aux;
    aux.tagi = Compression::multibyte_read(in);
    aux.tagj = Compression::multibyte_read(in);
    aux.tagk = Compression::multibyte_read(in);
    if(aux.tagk==999) aux.tagk=-999;
    forbid_rules.push_back(aux);
  }

  wcerr<<L"TaggerDataTrigram::read forbid_rules done\n";
  
  // array_tags
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    array_tags.push_back(Compression::wstring_read(in));
  }
  
  wcerr<<L"TaggerDataTrigram::read array_tags done\n";
  // tag_index
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    wstring tmp = Compression::wstring_read(in);    
    tag_index[tmp] = Compression::multibyte_read(in);
  }

  wcerr<<L"TaggerDataTrigram::read  tag_index done\n";
  // enforce_rules  
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    TEnforceAfterRule aux;
    aux.tagi = Compression::multibyte_read(in);
    for(int j = Compression::multibyte_read(in); j != 0; j--)
    {
      aux.tagsj.push_back(Compression::multibyte_read(in));
    }
    for(int k = Compression::multibyte_read(in); k != 0; k--)
    {
      aux.tagsk.push_back(Compression::multibyte_read(in));
    }
    enforce_rules.push_back(aux);
  }
  wcerr<<L"TaggerDataTrigram::read  enforce_rules done\n";

  // prefer_rules
  for(int i = Compression::multibyte_read(in); i != 0; i--)
  {
    prefer_rules.push_back(Compression::wstring_read(in));
  }

  wcerr<<L"TaggerDataTrigram::read  prefer_rules done\n";
  // constants
  constants.read(in);

  // output
  output.read(in); 

  // dimensions
  N = Compression::multibyte_read(in);
  M = Compression::multibyte_read(in);

  
  a = new double ** [N];
  b = new double ** [N];
  for(int i = 0; i != N; i++)
  {
    a[i] = new double *[N];
    b[i] = new double *[N];
    for(int j = 0; j != N; j++)
    {
    a[i][j] = new double [N];
    b[i][j] = new double [M];
    }
  }
   
  // read a
  for(int i = 0; i != N; i++)
  {
    for(int j = 0; j != N; j++)
    {
      for(int k = 0; k != N; k++){
        a[i][j][k] = EndianDoubleUtil::read(in);
        cerr<<"a["<<i<<"]["<<j<<"]["<<k<<"] = "<<a[i][j][k]<<"\n";
      }
    }
  }

  wcerr<<L"TaggerDataTrigram::read read a done\n";
  // initializing b matix
  for(int i = 0 ; i != N; i++)
  {
    for(int j = 0; j != N; j++)
    {
      for(int k = 0; k != M; k++)
      b[i][j][k] = 0.0;
    }
  }

  wcerr<<L"TaggerDataTrigram::read N="<<N<<L" M="<<M<<"\n";
  wcerr<<L"TaggerDataTrigram::read initialize b done\n";
  // read nonZERO values of b
  int nval = Compression::multibyte_read(in);

  for(; nval != 0; nval--)
  {
    int i = Compression::multibyte_read(in);
    int j = Compression::multibyte_read(in);
    int k = Compression::multibyte_read(in);
    cerr<<"b["<<i<<"]["<<j<<"]["<<k<<"] reading\n";
    b[i][j][k] = EndianDoubleUtil::read(in);
  }
  wcerr<<L"TaggerDataTrigram::read  values of b done\n";

  // read pattern list
  plist.read(in);
    
  // read discards on ambiguity
  discard.clear();

  int limit = Compression::multibyte_read(in);  
  if(feof(in))
  {
    return;
  }
  
  for(int i = 0; i < limit; i++)
  {
    discard.push_back(Compression::wstring_read(in));
  }
  wcerr<<L"TaggerDataTrigram::read over done\n";
}

void
TaggerDataTrigram::write(FILE *out)
{
  cerr << "TaggerDataTrigram::write" << endl;
  
  // open_class
  Compression::multibyte_write(open_class.size(), out);  
  int val = 0;
  for(set<TTag>::const_iterator it = open_class.begin(), limit = open_class.end();
      it != limit; it++)
  {
    Compression::multibyte_write(*it-val, out);    
    val = *it;
  }
  
  cerr << "TaggerDataTrigram::write open_class done" << endl;
  // forbid_rules
  Compression::multibyte_write(forbid_rules.size(), out);
  for(unsigned int i = 0, limit = forbid_rules.size(); i != limit; i++)
  {
    cerr << "TaggerDataTrigram:: inside forbid rules loop in\n" << endl;
    Compression::multibyte_write(forbid_rules[i].tagi, out);
    cerr << "TaggerDataTrigram:: inside forbid rules tagi doneloop\n" << endl;
    Compression::multibyte_write(forbid_rules[i].tagj, out);
    cerr << "TaggerDataTrigram:: inside forbid rules tagj done loop\n" << endl;
    wcerr<<L"tagk="<<forbid_rules[i].tagk<<endl;
    if(forbid_rules[i].tagk!=-999)
      Compression::multibyte_write(forbid_rules[i].tagk, out);
    else 
      Compression::multibyte_write(999, out);
    cerr << "TaggerDataTrigram:: inside forbid rules loop out\n" << endl;
  }
  
  cerr << "TaggerDataTrigram::write forbid_rules done" << endl;
  // array_tags
  Compression::multibyte_write(array_tags.size(), out);
  for(unsigned int i = 0, limit = array_tags.size(); i != limit; i++)
  {
    Compression::wstring_write(array_tags[i], out);
  }

  cerr << "TaggerDataTrigram::write array_tags done" << endl;
  // tag_index
  Compression::multibyte_write(tag_index.size(), out);
  for(map<wstring, int, Ltstr>::iterator it = tag_index.begin(), limit = tag_index.end();
      it != limit; it++)
  {
    Compression::wstring_write(it->first, out);
    Compression::multibyte_write(it->second, out);
  }
  
  cerr << "TaggerDataTrigram::write tag_index done" << endl;
  // enforce_rules
  Compression::multibyte_write(enforce_rules.size(), out);
  for(unsigned int i = 0, limit = enforce_rules.size(); i != limit; i++)
  {
    Compression::multibyte_write(enforce_rules[i].tagi, out);
    Compression::multibyte_write(enforce_rules[i].tagsj.size(), out);
    for(unsigned int j = 0, limit2 = enforce_rules[i].tagsj.size(); j != limit2; j++)
    {
      Compression::multibyte_write(enforce_rules[i].tagsj[j], out);
    }
    Compression::multibyte_write(enforce_rules[i].tagsk.size(), out);
    for(unsigned int k = 0, limit2 = enforce_rules[i].tagsk.size(); k != limit2; k++)
    {
      Compression::multibyte_write(enforce_rules[i].tagsk[k], out);
    }
  }

  cerr << "TaggerDataTrigram::write  enforce_rules done" << endl;
  // prefer_rules
  Compression::multibyte_write(prefer_rules.size(), out);
  for(unsigned int i = 0, limit = prefer_rules.size(); i != limit; i++)
  {
    Compression::wstring_write(prefer_rules[i], out);
  }
  
  cerr << "TaggerDataTrigram::write prefer_rules done" << endl;
  // constants
  constants.write(out);  

  // output
  output.write(out);

  // a matrix
  Compression::multibyte_write(N, out);
  Compression::multibyte_write(M, out);
  for(int i = 0; i != N; i++)
  {
    for(int j = 0; j != N; j++)
    {
      for(int k = 0; k != N; k++)
      EndianDoubleUtil::write(out, a[i][j][k]);
    }
  }
  cerr << "TaggerDataTrigram::write a matrix done" << endl;

  // b matrix, writing only useful values
  
  int nval = 0;
  for(int i = 0; i != N; i++)
  {
    for(int j = 0; j != N; j++)
    {
      for(int k = 0; k != M; k++)
      {
        if(output[k].find(j) != output[k].end())
        {
	  nval++;
        }
      }
    }
  }

  Compression::multibyte_write(nval, out);
  for(int i = 0; i != N; i++)
  {
    for(int j = 0; j != N; j++)
    {
      for(int k = 0; k != M; k++)
      {
        if(output[k].find(j) != output[k].end())
        {
	  Compression::multibyte_write(i, out);
	  Compression::multibyte_write(j, out);
	  Compression::multibyte_write(k, out);
	  EndianDoubleUtil::write(out, b[i][j][k]);
        }
      }
    }
  }  
  cerr << "TaggerDataTrigram::write  b matrix done" << endl;
  
  // write pattern list
  plist.write(out);
  
  // write discard list
  
  if(discard.size() != 0)
  {
    Compression::multibyte_write(discard.size(), out);
    for(unsigned int i = 0, limit = discard.size(); i != limit; i++)
    {
      Compression::wstring_write(discard[i], out);
    }
  }  
  cerr << "TaggerDataTrigram::write discard list done" << endl;
}

void
TaggerDataTrigram::addDiscard(wstring const &tags)
{
  discard.push_back(tags);
}