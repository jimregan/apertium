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
/**
 *  Light Sliding-Window Part of Speech Tagger (LSWPoST) implementation (source)
 *
 *  @author   Gang Chen - pkuchengang@gmail.com
 */


#include <apertium/lswpost.h>
#include <apertium/tagger_utils.h>
#include  "apertium_config.h"
#include <apertium/unlocked_cstdio.h>
#include <lttoolbox/compression.h>

#ifdef WIN32
#define isnan(n) _isnan(n)
#define isinf(n) (!_finite(n))
#endif

#ifdef __clang__
#undef __GNUC__
#endif

#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <apertium/string_utils.h>

using namespace std;
using namespace Apertium;
using namespace tagger_utils;

LSWPoST::LSWPoST(TaggerDataLSW *t) {
  this->tdlsw = t;
  debug=false;
  show_sf=false;
  null_flush = false;
  eos = (tdlsw->getTagIndex())[L"TAG_SENT"];  
}

LSWPoST::~LSWPoST() {
}

void
LSWPoST::set_eos(TTag t) { 
  eos = t; 
} 

void
LSWPoST::set_debug(bool d) { 
  debug = d; 
} 

void
LSWPoST::set_show_sf(bool sf) { 
  show_sf = sf; 
}
 
void
LSWPoST::setNullFlush(bool nf) {
  null_flush = nf;
}

void
LSWPoST::init_probabilities(FILE *ftxt) {

  int N = tdlsw->getN();
  int nw = 0;

  vector<vector<vector<double> > > para_matrix(N, vector<vector<double> >(N, vector<double>(N, 0)));

  Collection &output = tdlsw->getOutput();

  MorphoStream morpho_stream(ftxt, true, tdlsw);
  
  set<TTag> tags_left, tags, tags_right;
  
  TaggerWord *word_left = new TaggerWord(); // word_left
  word_left->add_tag(eos, L"sent", tdlsw->getPreferRules());
  TaggerWord *word = morpho_stream.get_next_word(); // word
  if (morpho_stream.getEndOfFile()) {
    delete word_left;
    delete word;
    return;
  }
  TaggerWord *word_right = morpho_stream.get_next_word(); // word_right

  //We count each element of the para matrix
  while (word_right != NULL) {
    if (++nw % 10000 == 0) {
      wcerr << L'.' << flush;
    }

    tags_left = word_left->get_tags();
    if (tags_left.size() == 0) { //This is an unknown word
      tags_left = tdlsw->getOpenClass();
    }
    tags = word->get_tags();
    if (tags.size() == 0) {
      tags = tdlsw->getOpenClass();
    }
    tags_right = word_right->get_tags();
    if (tags_right.size() == 0) {
      tags_right = tdlsw->getOpenClass();
    }

    if (output.has_not(tags_left) || output.has_not(tags) || output.has_not(tags_right)) {
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors += L"Word '" + word->get_superficial_form() + L"' not found in the dictionary.\n";
      errors += L"New ambiguity class: " + word->get_string_tags() + L"\n";
      errors += L"Take a look at the dictionary and at the training corpus. Then, retrain.";
      fatal_error(errors);
    }
    set<TTag>::iterator iter, iter_left, iter_right;

    for (iter = tags.begin(); iter != tags.end(); ++iter) {
      for (iter_left = tags_left.begin(); iter_left != tags_left.end(); ++iter_left) {
        for (iter_right = tags_right.begin(); iter_right != tags_right.end(); ++iter_right) {
          para_matrix[*iter_left][*iter][*iter_right] +=
              1.0 / (tags_left.size() * tags.size() * tags_right.size());
        }
      }
    }

    delete word_left;
    word_left = word;
    word = word_right;
    word_right = morpho_stream.get_next_word();
  }
  delete word_left;
  delete word;

  //tdlsw->setLSWPoSTProbabilities(N, M, para_matrix);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        tdlsw->getD()[i][j][k] = para_matrix[i][j][k];
      }
    }
  }

  wcerr << L"\n";
}

void
LSWPoST::apply_rules() {

  vector<TForbidRule> &forbid_rules = tdlsw->getForbidRules();
  vector<TEnforceAfterRule> &enforce_rules = tdlsw->getEnforceRules();
  int N = tdlsw->getN();

  /** Forbid Rules
   *  for any sequence that contains a forbid rule,
   *  the prob is set to ZERO.
   */
  for (size_t r = 0; r < forbid_rules.size(); ++r) {
    TTag tagi = forbid_rules[r].tagi;
    TTag tagj = forbid_rules[r].tagj;
    for (int k = 0; k < N; ++k) {
      //         left mid right
      tdlsw->getD()[tagi][tagj][k] = 0;
      tdlsw->getD()[k][tagi][tagj] = 0;
    }
  }


  /** Enforce Rules
   *  for any sequence that doesn't satisfy a enforce rule,
   *  the prob is set to ZERO.
   */
  for (size_t r = 0; r < enforce_rules.size(); ++r) {
    TTag tagi = enforce_rules[r].tagi;
    vector<TTag> tagsj = enforce_rules[r].tagsj;
    for (TTag tagj = 0; tagj < N; ++tagj) {
      bool found = false;
      for (size_t j = 0; j < tagsj.size(); ++j) {
        if (tagj == tagsj[j]) {
          found = true;
          break;
        }
      }
      if (!found) {
        for (int k = 0; k < N; ++k) {
          //         left mid right
          tdlsw->getD()[tagi][tagj][k] = 0;
          tdlsw->getD()[k][tagi][tagj] = 0;
        }
      }
    } // for tagj
  } // for r
}

void
LSWPoST::read_dictionary(FILE *fdic) {
  int i, k, nw = 0;
  TaggerWord *word = NULL;
  set<TTag> tags;
  Collection &output = tdlsw->getOutput();

  MorphoStream morpho_stream(fdic, true, tdlsw);

  // In the input dictionary there must be all punctuation marks, including the end-of-sentece mark

  word = morpho_stream.get_next_word();

  while (word) {
    if (++nw % 10000 == 0)
      wcerr << L'.' << flush;

    tags = word->get_tags();

    if (tags.size() > 0)
      k = output[tags];

    delete word;
    word = morpho_stream.get_next_word();
  }
  wcerr << L"\n";

  // OPEN AMBIGUITY CLASS
  // It contains all tags that are not closed.
  // Unknown words are assigned the open ambiguity class
  k = output[tdlsw->getOpenClass()];

  int N = (tdlsw->getTagIndex()).size();

  // Create ambiguity class holding one single tag for each tag.
  // If not created yet
  for (i = 0; i != N; i++) {
    set < TTag > amb_class;
    amb_class.insert(i);
    k = output[amb_class];
  }

  wcerr << N << L" states\n";
  tdlsw->setProbabilities(N);
}

void
LSWPoST::train(FILE *ftxt) {
  int N = tdlsw->getN();
  int nw = 0;

  vector<vector<vector<double> > > para_matrix_new(N, vector<vector<double> >(N, vector<double>(N, 0)));

  set<TTag> tags_left, tags, tags_right;

  Collection &output = tdlsw->getOutput();

  MorphoStream morpho_stream(ftxt, true, tdlsw);

  TaggerWord *word_left = new TaggerWord(); // word_left
  word_left->add_tag(eos, L"sent", tdlsw->getPreferRules());
  TaggerWord *word = morpho_stream.get_next_word(); // word
  if (morpho_stream.getEndOfFile()) {
    delete word_left;
    delete word;
    return;
  }
  TaggerWord *word_right = morpho_stream.get_next_word(); // word_right

  while (word_right) {

    //wcerr<<L"Enter para continuar\n";
    //getchar();
    if (++nw % 10000 == 0) {
      wcerr << L'.' << flush;
    }

    tags_left = word_left->get_tags();
    if (tags_left.size() == 0) { //This is an unknown word
      tags_left = tdlsw->getOpenClass();
    }
    tags = word->get_tags();
    if (tags.size() == 0) { //This is an unknown word
      tags = tdlsw->getOpenClass();
    }
    tags_right = word_right->get_tags();
    if (tags_right.size() == 0) { //This is an unknown word
      tags_right = tdlsw->getOpenClass();
    }

    if (output.has_not(tags_left) || output.has_not(tags) || output.has_not(tags_right)) {
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors += L"Word '" + word->get_superficial_form() + L"' not found in the dictionary.\n";
      errors += L"New ambiguity class: " + word->get_string_tags() + L"\n";
      errors += L"Take a look at the dictionary and at the training corpus. Then, retrain.";
      fatal_error(errors);
    }

    double normalization = 0;
    set<TTag>::iterator iter, iter_left, iter_right;

    for (iter = tags.begin(); iter != tags.end(); ++iter) {
      for (iter_left = tags_left.begin(); iter_left != tags_left.end(); ++iter_left) {
        for (iter_right = tags_right.begin(); iter_right != tags_right.end(); ++iter_right) {
          normalization += tdlsw->getD()[*iter_left][*iter][*iter_right];
        }
      }
    }

    for (iter = tags.begin(); iter != tags.end(); ++iter) {
      for (iter_left = tags_left.begin(); iter_left != tags_left.end(); ++iter_left) {
        for (iter_right = tags_right.begin(); iter_right != tags_right.end(); ++iter_right) {
          if (normalization > ZERO) {
            para_matrix_new[*iter_left][*iter][*iter_right] += 1.0 / normalization;
          }
        }
      }
    }

    delete word_left;
    word_left = word;
    word = word_right;
    word_right = morpho_stream.get_next_word();
  }
  delete word_left;
  delete word;

  //tdlsw-setLSWPoSTProbabilities(N, M, (double ***)para_matrix_new);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        tdlsw->getD()[i][j][k] = tdlsw->getD()[i][j][k] * para_matrix_new[i][j][k];
      }
    }
  }
}

void
LSWPoST::print_para_matrix() {
  wcout << L"para matrix D\n----------------------------\n";
  for (int i = 0; i < tdlsw->getN(); ++i) {
    for (int j = 0; j < tdlsw->getN(); ++j) {
      for (int k = 0; k < tdlsw->getN(); ++k) {
        wcout << L"D[" << i << L"][" << j << L"][" << k << L"] = "
            << tdlsw->getD()[i][j][k] << "\n";
      }
    }
  }
}

void 
LSWPoST::tagger(FILE *in, FILE *out, bool show_all_good_first) {
  set<TTag> tags_left, tags, tags_right;

  MorphoStream morpho_stream(in, debug, tdlsw);
  morpho_stream.setNullFlush(null_flush);                      
  
  Collection &output = tdlsw->getOutput();
  
  TaggerWord *word_left = new TaggerWord(); // word_left
  word_left->add_tag(eos, L"sent", tdlsw->getPreferRules());
  word_left->set_show_sf(show_sf);
  TaggerWord *word = morpho_stream.get_next_word(); // word
  if (morpho_stream.getEndOfFile()) {
    delete word_left;
    delete word;
    return;
  }
  word->set_show_sf(show_sf);
  TaggerWord *word_right = morpho_stream.get_next_word(); // word_right
  word_right->set_show_sf(show_sf);

  wstring micad;

  while (word_right) {
    tags_left = word_left->get_tags();
    if (tags_left.size() == 0) {
      tags_left = tdlsw->getOpenClass();
    }
    tags = word->get_tags();
    if (tags.size() == 0) {
      tags = tdlsw->getOpenClass();
    }
    tags_right = word_right->get_tags();
    if (tags_right.size() == 0) {
      tags_right = tdlsw->getOpenClass();
    }
    if (output.has_not(tags_right)) {
      if (debug) {
        wstring errors;
        errors = L"A new ambiguity class was found. \n";
        errors+= L"Retraining the tagger is necessary so as to take it into account.\n";
        errors+= L"Word '"+word_right->get_superficial_form()+L"'.\n";
        errors+= L"New ambiguity class: "+word_right->get_string_tags()+L"\n";
        fatal_error(errors);
      }
      tags_right = find_similar_ambiguity_class(tags_right);
    }

    double max = -1;
    TTag tag_max = *tags.begin();
    set<TTag>::iterator iter, iter_left, iter_right;
    for (iter = tags.begin(); iter != tags.end(); ++iter) {
      double n = 0;
      for (iter_left = tags_left.begin(); iter_left != tags_left.end(); ++iter_left) {
        for (iter_right = tags_right.begin(); iter_right != tags_right.end(); ++iter_right) {
          n += tdlsw->getD()[*iter_left][*iter][*iter_right];
        }
      }
      if (n > max) {
        max = n;
        tag_max = *iter;
      }
    }

    micad = word->get_lexical_form(tag_max, (tdlsw->getTagIndex())[L"TAG_kEOF"]);
    fputws_unlocked(micad.c_str(), out);
    if (morpho_stream.getEndOfFile()) {
      if (null_flush) {
        fputwc_unlocked(L'\0', out);
      }
      fflush(out);
      morpho_stream.setEndOfFile(false);
    }
  
    delete word_left;
    word_left = word;
    word = word_right;
    word_right = morpho_stream.get_next_word();
    if (word_right != NULL) {
      word_right->set_show_sf(show_sf);
    }
  }
  delete word_left;
  delete word;
}

set<TTag>                                                                                     
LSWPoST::find_similar_ambiguity_class(set<TTag> c) {
  int size_ret = -1;                                                                          
  set<TTag> ret=tdlsw->getOpenClass(); // return open-class as default, if no better is found.
  bool skip_class;                                                                           
  Collection &output = tdlsw->getOutput();                                                       
                                                                                              
  for(int k=0; k<output.size(); k++) {                                                           
    if ((((int)output[k].size())>((int)size_ret)) && (((int)output[k].size())<((int)c.size()))) {    
      skip_class=false;                                                                      
      // Test if output[k] is a subset of class                                               
      for(set<TTag>::const_iterator it=output[k].begin(); it!=output[k].end(); it++) {        
        if (c.find(*it)==c.end()) { 
           skip_class=true; //output[k] is not a subset of class                             
           break;
        }  
      } 
      if (!skip_class) {                                                                     
        size_ret = output[k].size();                                                          
             ret = output[k];
      }      
    } 
  } 
  return ret;                                                                                 
} 

