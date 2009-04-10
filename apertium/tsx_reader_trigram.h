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
#ifndef _TSXREADERTRIGRAM_
#define _TSXREADERTRIGRAM_

#include <apertium/constant_manager.h>
#include "tagger_data_trigram.h"
#include <apertium/ttag.h>
#include <lttoolbox/pattern_list.h>
#include <lttoolbox/ltstr.h>

#include <libxml/xmlreader.h>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

class TSXReaderTrigram
{
private:
  xmlTextReaderPtr reader;  
  set<TTag> *open_class;
  vector<TForbidRule> *forbid_rules;
  map<wstring, TTag, Ltstr> *tag_index;
  vector<wstring> *array_tags;
  vector<TEnforceAfterRule> *enforce_rules;
  vector<wstring> *prefer_rules;
  PatternList *plist;
  ConstantManager *constants;
  TaggerDataTrigram tdata;

  int type;
  wstring name;

  wstring attrib(wstring const &name);

  void parseError(wstring const &message);
  void newTagIndex(wstring const &tag);
  void newDefTag(wstring const &tag);
  void newConstant(wstring const &constant);
  void procDefLabel();
  void procDefMult();
  void procDiscardOnAmbiguity();
  void procTagset();
  void procForbid();
  void procLabelSequence();
  void procEnforce();
  void procPreferences();
  void copy(TSXReaderTrigram const &o);
  void destroy();
  void clearTagIndex();

  void step();
public:
  TSXReaderTrigram();
  ~TSXReaderTrigram();
  TSXReaderTrigram(TSXReaderTrigram const &o);
  TSXReaderTrigram & operator =(TSXReaderTrigram const &o);

  void read(string const &filename);
  void write(string const &filename);
  TaggerDataTrigram & getTaggerData();
};

#endif
