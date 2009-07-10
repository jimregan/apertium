
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
/*
 *  2nd order hidden Markov model (HMM) implementation (source)
 *
 *  @author	Felipe Sánchez-Martínez - fsanchez@dlsi.ua.es
 */

//#include <apertium/hmm2.h>
#include "hmm2.h"
#include <apertium/tagger_utils.h>
#include <apertium/unlocked_cstdio.h>
#include <lttoolbox/compression.h>
#include "smooth_utils_trigram.h"

#ifdef WIN32
#define isnan(n) _isnan(n)
#define isinf(n) (!_finite(n))
#endif

#include <stdio.h>
#include <limits>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <apertium/string_utils.h>

using namespace Apertium;
using namespace tagger_utils;

HMM2::HMM2(TaggerDataTrigram *t)
{
  //cerr<<"HMM2::constructor\n";
  this->td = t;

  debug=false;
  show_sf=false;
  eos = (td->getTagIndex())[L"TAG_SENT"];  
}

HMM2::~HMM2()
{
}

void
HMM2::init()
{
}

void
HMM2::set_eos(TTag t) 
{ 
  eos = t; 
} 

void
HMM2::set_debug(bool d)
{ 
  debug = d; 
} 

void
HMM2::set_show_sf(bool sf)
{ 
  show_sf = sf; 
} 


void 
HMM2::read_ambiguity_classes(FILE *in) 
{
  while(in)
  {
    int ntags = Compression::multibyte_read(in);

    if(feof(in))
    {
      break;
    }
    set<TTag> ambiguity_class;

    for(; ntags != 0; ntags--)
    {
      ambiguity_class.insert(Compression::multibyte_read(in));
    }
    
    if(ambiguity_class.size() != 0)
    {
      td->getOutput().add(ambiguity_class);
    }     
  }
  
  td->setProbabilities(td->getTagIndex().size(), td->getOutput().size());
}

void 
HMM2::write_ambiguity_classes(FILE *out) 
{
  for(int i=0, limit = td->getOutput().size(); i != limit; i++) 
  {
    set<TTag> const &ac = (td->getOutput())[i];
    Compression::multibyte_write(ac.size(), out);
    for(set<TTag>::const_iterator it = ac.begin(), limit2 = ac.end();
        it != limit2; it++)
    {
      Compression::multibyte_write(*it, out);
    }
  } 
}  

void 
HMM2::read_probabilities(FILE *in)
{
  td->read(in);
}

void 
HMM2::write_probabilities(FILE *out)
{
  td->write(out);  
}  

void 
HMM2::init_probabilities_kupiec (FILE *is, int corpus_length, string savecountsfile)
{
  //cerr<<"HMM2::Kupiec\n";
  int N = td->getN();
  int M = td->getM();
  int i, j, k, k1, k2, k3, nw=0;
  map<int, double> classes_ocurrences; // M
  map<int, map<int, double> > classes_pair_ocurrences; // MxM
  map<int, map<int,  map<int, double> > > classes_triple_ocurrences; // MxMxM

  map<int, double> tags_count; // N
  map<int, double> tags_count_for_emis; // N
  map<int,  map<int, double> > tags_pair_for_emis2; // NxN
  map<int, map<int, double> > tags_pairs; //NxN
  map<int, map<int,  map<int, double> > > tags_triple; //NxNxN
  map<int, map<int, double> > emis; //NxM
  map<int, map<int,  map<int, double> > > emis2; //NxNxM


  Collection &output = td->getOutput();
 
  MorphoStreamTrigram lexmorfo(is, true, td);
  
  TaggerWord *word=NULL;
/*
  for(k=0; k<M; k++) {
    classes_ocurrences[k]=0; 
    for (k2=0; k2<M; k2++){
      classes_pair_ocurrences[k][k2]=0;
        for (k3=0; k3<M; k3++)
	  classes_triple_ocurrences[k][k2][k3]=0;
      }
  }
*/
  set<TTag> tags;
  tags.clear();
  tags.insert(eos);  
  k1=output[tags]; //The first tag (ambiguity class) seen is the end-of-sentence
  classes_ocurrences[k1]++;
  
  //We count for each ambiguity class the number of ocurrences
  word = lexmorfo.get_next_word();

    if (++nw%10000==0) wcerr<<L'.'<<flush;
   
    tags=word->get_tags();

    if (tags.size()==0) { //This is an unknown word
      tags = td->getOpenClass();
    }
    else if (output.has_not(tags)) {
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors+= L"Word '"+word->get_superficial_form()+L"' not found in the dictionary.\n";
      errors+= L"New ambiguity class: "+word->get_string_tags()+L"\n";
      errors+= L"Take a look at the dictionary and at the training corpus. Then, retrain.";
      fatal_error(errors);
    }

    k2=output[tags];

    classes_ocurrences[k1]++;
    classes_pair_ocurrences[k1][k2]++;  //k1 followed by k2
    delete word;
    word=lexmorfo.get_next_word();

  while((word)) {
    //wcerr<<"Kupiec word: "<<word->get_superficial_form()<<"\n";
    if (++nw%10000==0) wcerr<<L'.'<<flush; 
    
    tags=word->get_tags();

    if (tags.size()==0) { //This is an unknown word
      tags = td->getOpenClass();
    }
    else if (output.has_not(tags)) { 
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors+= L"Word '"+word->get_superficial_form()+L"' not found in the dictionary.\n";
      errors+= L"New ambiguity class: "+word->get_string_tags()+L"\n";
      errors+= L"Take a look at the dictionary and at the training corpus. Then, retrain.";      
      fatal_error(errors);      
    }    

    k3=output[tags];

    classes_ocurrences[k1]++;
    classes_pair_ocurrences[k1][k2]++;  //k1 followed by k2
    classes_triple_ocurrences[k1][k2][k3]++;  //k1 followed by k2 followed by k3
    delete word;
    word=lexmorfo.get_next_word();

    k1=k2;
    k2=k3;

    if((corpus_length>0)&&(nw>=corpus_length)&&(tags.size()==1)&&((*(tags.begin()))==eos)) {
      //En of training
      cerr<<nw<<" ";
      break;
    }
  }
 
  
  //wcerr<<"Kupiec word: Out of while loop\n";

  //Estimation of the number of time each tags occurs in the training text
  for(i=0; i<N; i++) {  
    tags_count_for_emis[i]=0;
    for(k=0; k<M;  k++) { 
        if(output[k].find(i) != output[k].end())
        tags_count_for_emis[i] += classes_ocurrences[k]/((double)output[k].size());
    }
  }

 
  //wcerr<<"Kupiec: tags_count_for_emis done\n";
  //Estimation of the number of times each tag pair occurs
/*
  for(i=0; i<N; i++){
    tags_count[i]=0;
    for(j=0; j<N; j++){
      //tags_pair_estimate[i][j]=0;
      tags_pairs[i][j]=0;
      tags_pair_for_emis2[i][j]=0;
      for(k=0;k<N;k++){
        tags_triple[i][j][k]=0;
      }        
      for(k=0;k<M;k++){ 
        emis2[i][j][k]=0;
      }        
    }
  }
*/
  //wcerr<<"Kupiec: tags_count etc zeroed done\n";
  set<TTag> tags1, tags2, tags3;
  set<TTag>::iterator itag1, itag2, itag3;
  for(k1=0; k1<M; k1++) {
    tags1=output[k1];
    for(k2=0; k2<M; k2++) {
      tags2=output[k2];
      for(k3=0; k3<M; k3++) {
        tags3=output[k3];
        double nocurrences=classes_triple_ocurrences[k1][k2][k3]/((double)(tags1.size()*tags2.size()*tags3.size()));
        for (itag1=tags1.begin(); itag1!=tags1.end(); itag1++) {
	  if (((*itag1)<0)||((*itag1)>=N))
            cerr<<"Error: Tag "<<*itag1<<" out of range\n";
          for (itag2=tags2.begin(); itag2!=tags2.end(); itag2++) {
            if (((*itag2)<0)||((*itag2)>=N))
              cerr<<"Error: Tag "<<*itag2<<" out of range\n";
            for (itag3=tags3.begin(); itag3!=tags3.end(); itag3++) {
              if (((*itag3)<0)||((*itag3)>=N))
                cerr<<"Error: Tag "<<*itag3<<" out of range\n";
	      tags_triple[*itag1][*itag2][*itag3] += nocurrences;
              tags_pairs[*itag1][*itag2]+=nocurrences;
              tags_count[*itag1]+=nocurrences;
            }
	  }
        }
      }
    }
  }

  //wcerr<<"Kupiec: tags_count etc done done\n";

  //Estimation of the number of times each ambiguity class is emitted
  //from each tag
  for(i=0; i<N; i++) {
    for(k=0; k<M; k++)  {
      if (output[k].find(i)!=output[k].end()) {
        emis[i][k]=classes_ocurrences[k]/((double)output[k].size());
      }
    }
  }
  //wcerr<<"Kupiec: emis done\n";

  //Estimation of the number of times each ambiguity class is emitted
  //from each tag pair
  

  for(k1=0; k1<M; k1++) {
    tags1=output[k1];
    for(k2=0; k2<M; k2++) {
      tags2=output[k2];
      double nocurrences=classes_pair_ocurrences[k1][k2]/((double)(tags1.size()*tags2.size()));
      for (itag1=tags1.begin(); itag1!=tags1.end(); itag1++) {
        if (((*itag1)<0)||((*itag1)>=N))
          cerr<<"Error: Tag "<<*itag1<<" out of range\n";
        for (itag2=tags2.begin(); itag2!=tags2.end(); itag2++) {
          if (((*itag2)<0)||((*itag2)>=N))
            cerr<<"Error: Tag "<<*itag2<<" out of range\n";
          emis2[*itag1][*itag2][k2] += nocurrences;
          tags_pair_for_emis2[*itag1][*itag2] += nocurrences;
        }
      }
    }
  }
  //wcerr<<"Kupiec: emis2  done\n";
/*
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      for(k=0; k<M; k++)  {
        if (output[k].find(j)!=output[k].end()) {
          for(k1=0;k1<M;k1++){
            if (output[k1].find(i)!=output[k1].end()) 
              emis2[i][j][k]+=classes_pair_ocurrences[k1][k]/((double)(output[k1].size()*output[k].size()));
              cerr<<"emis2["<<i<<"]["<<j<<"]["<<k<<"] ="<<emis2[i][j][k]<<"\n";
          }
        }
        tags_pair_for_emis2[i][j] += emis2[i][j][k];
      }
    }
  }
*/

//jun23
/*  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      for(k=0; k<M; k++)  {
        tags_pair_for_emis2[i][j] += emis2[i][j][k];
      }
    }
  }
*/
  //wcerr<<"Kupiec: tags_pair_for_emis2 done\n";


  if (savecountsfile!="") {
    cerr<<"Saving counts to file '"<<savecountsfile<<"' ... "<<flush;
    SmoothUtilsTrigram::save_counts(*td, savecountsfile, tags_count, tags_pairs, tags_triple, classes_ocurrences, emis2, emis, tags_pair_for_emis2, tags_count_for_emis);
    cerr<<"done.\n";
  }

  SmoothUtilsTrigram::calculate_smoothed_parameters(*td, tags_count, tags_pairs, tags_triple, classes_ocurrences, emis2, emis, tags_pair_for_emis2, tags_count_for_emis, nw);
  //SmoothUtilsTrigram::calculate_smoothed_parameters(*td, tags_count_for_emis, tags_pair_for_emis2, tags_triple, classes_ocurrences, emis2, emis, tags_pair_for_emis2, tags_count_for_emis, nw);
  cerr<<" done\n";
}

void 
HMM2::init_probabilities_from_tagged_text(FILE *ftagged, FILE *funtagged, string savecountsfile) {
  //cerr<<"HMM2::train_tagged\n";
  int i, j, k, nw=0;
  int N = td->getN();
  int M = td->getM();

  map<int, map<int, double> > tags_pair; //NxN
  map<int, map<int, map<int, double> > > tags_triple; //NxNxN
  map<int, map<int, double> > emis;     //NxM
  map<int, map<int, map<int, double> > > emis2;     //NxNxM
  map<int, double> ambclass_count;      //M
  map<int, double> tags_count;          //N
  map<int, double> tags_count_for_emis; //N
  map<int, map<int, double> > tags_pair_for_emis2; //NxN
  
  MorphoStreamTrigram stream_tagged(ftagged, true, td);
  MorphoStreamTrigram stream_untagged(funtagged, true, td);
  
  TaggerWord *word_tagged=NULL, *word_untagged=NULL;
  Collection &output = td->getOutput();

  
  set<TTag> tags;
 
  // Init counters - each event appears at least once. 
  // Espected likelihood estimate (ELE) with a fixed initial count of 1
/*  for(i=0; i<N; i++) {
    tags_count[i]=0;
    tags_count_for_emis[i]=0;

    for(j=0; j<N; j++){
      tags_pair[i][j]=0;
      tags_pair_for_emis2[i][j]=0;

      for(k=0;k<N;k++)
	tags_triple[i][j][k]=0;	  
    }
  }

  for(k=0; k<M; k++) {
    ambclass_count[k]=0;
    for(i=0; i<N; i++) {
      if (output[k].find(i)!=output[k].end())
        emis[i][k] = 0;
      for(j=0; j<N; j++) {
        if (output[k].find(j)!=output[k].end())
	emis2[i][j][k]=0;
      }
    }  
  }
*/ 
  TTag tag1, tag2, tag3;  
  tag1 = eos; // The first seen tag is the end-of-sentence tag
  tag2 = eos;
  tag3 = eos; //Doesn't really matter, but I need it initialized
  
  word_tagged = stream_tagged.get_next_word();
  word_untagged = stream_untagged.get_next_word();
  while(word_tagged) {
   // wcerr<<*word_tagged;
   // wcerr<<L" -- "<<*word_untagged<<L"\n"; 

    if (word_tagged->get_superficial_form()!=word_untagged->get_superficial_form()) {              
      wcerr<<L"\nTagged text (.tagged) and analyzed text (.untagged) streams are not aligned.\n";
      wcerr<<L"Take a look at tagged text (.tagged).\n";
      wcerr<<L"Perhaps this is caused by a multiword unit that is not a multiword unit in one of the two files.\n";
      wcerr<<*word_tagged<<L" -- "<<*word_untagged<<L"\n"; 
      exit(1);
    }

    if (++nw%100==0) wcerr<<L'.'<<flush; 
    
    tag3 = tag2;
    tag2 = tag1;
   
    if (word_untagged==NULL) {
      wcerr<<L"word_untagged==NULL\n";
      exit(1);
    }

    if (word_tagged->get_tags().size()==0) // Unknown word
    {
	tag1 = -1; 
      //tag1=td->getTagIndex()[L"TAG_kUNDEF"];
    }
    else if (word_tagged->get_tags().size()>1) // Ambiguous word
      wcerr<<L"Error in tagged text. An ambiguous word was found: "<<word_tagged->get_superficial_form()<<L"\n";
    else
      tag1 = *(word_tagged->get_tags()).begin();


    if ((tag1>=0) && (tag2>=0) && (tag3>0)){
      tags_triple[tag3][tag2][tag1]++;
      tags_pair[tag3][tag2]++;
      tags_count[tag3]++;
    }
    

    if (word_untagged->get_tags().size()==0) { // Unknown word
      tags = td->getOpenClass();
    }
    else if (output.has_not(word_untagged->get_tags())) { //We are training, there is no problem
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors+= L"Word '"+word_untagged->get_superficial_form()+L"' not found in the dictionary.\n";
      errors+= L"New ambiguity class: "+word_untagged->get_string_tags()+L"\n";
      errors+= L"Take a look at the dictionary, then retrain.";
      fatal_error(errors);      
    }    
    else {
      tags = word_untagged->get_tags();
    }

    //if(tag1>=N) cerr<<"Whoops, tag greater than "<<N<<": "<<tag1<<" k="<<k<<"\n";
    //if(tag1<0) cerr<<"Whoops, tag less than 0: "<<tag1<<" k="<<k<<"\n";
    k=output[tags];
    if(tag2>=0 && tag1>=0){
      if (output[k].find(tag1)!=output[k].end()) {
	emis2[tag2][tag1][k]++;
        emis[tag1][k]++;
        tags_count_for_emis[tag1]++;
        tags_pair_for_emis2[tag2][tag1]++;
        ambclass_count[k]++;
	//cerr<<"Warning: updating amb counts\n";
      } else {
        cerr<<"Warning: Ambiguity class "<<k<<" is emmited from tag "<<tag1<<" but it should not\n";
      }
    }
                   
    delete word_tagged;
    word_tagged=stream_tagged.get_next_word();
    delete word_untagged;
    word_untagged=stream_untagged.get_next_word();       
  }
  //tags_pair[tag2][tag1]++;
  //tags_count[tag2]++;
  //tags_count[tag1]++;

  SmoothUtilsTrigram::calculate_smoothed_parameters(*td, tags_count, tags_pair ,tags_triple,  ambclass_count, emis2, emis, tags_pair_for_emis2, tags_count_for_emis, nw);

  if (savecountsfile!="") {
    cerr<<"Saving counts to file '"<<savecountsfile<<"' ... "<<flush;
    SmoothUtilsTrigram::save_counts(*td, savecountsfile, tags_count, tags_pair, tags_triple, ambclass_count, emis2, emis, tags_pair_for_emis2,tags_count_for_emis);
    cerr<<"done.\n";
  }

  cerr<<"Number of words processed: "<<nw<<"\n";
}
  
void
HMM2::apply_rules()
{

  //cerr<<"HMM2::applyrules\n";
  vector<TForbidRule> &forbid_rules = td->getForbidRules();
  vector<TEnforceAfterRule> &enforce_rules = td->getEnforceRules();
  int N = td->getN();
  int i, j, j2, k, k2;
  bool found;
 
  for(i=0; i<(int) forbid_rules.size(); i++) {
    if(forbid_rules[i].tagk==-999){
    //cerr<<"HMM2::applyrules two label-item forbid rule\n";
   
      //if only two label-items specified in forbid rule
      for(k=0; k<N; k++){ 
        (td->getA())[forbid_rules[i].tagi][forbid_rules[i].tagj][k] = ZERO;
        (td->getA())[k][forbid_rules[i].tagi][forbid_rules[i].tagj] = ZERO;
      }
    }
    else
      //if three label-items specified in forbid rule
      (td->getA())[forbid_rules[i].tagi][forbid_rules[i].tagj][forbid_rules[i].tagk] = ZERO;
    //cerr<<"HMM2::applyrules 3 label-item forbid rule\n";
  }
  //cerr<<"HMM2::applyrules forbid rules done\n";

  for(i=0; i<(int) enforce_rules.size(); i++) {
    for(j=0; j<N; j++) {
      for(k=0; k<N; k++) {
        found = false;
        for (j2=0; j2<(int) enforce_rules[i].tagsj.size(); j2++) {
	  if (enforce_rules[i].tagsj[j2]==j) {              //enforce rules with tagsj only
            if(enforce_rules[i].tagsk.size()==0){
	      found = true;
	      (td->getA())[k][enforce_rules[i].tagi][j] = ZERO; //any tag, tagi, tagsj
            }
	    for (k2=0; k2<(int) enforce_rules[i].tagsk.size(); k2++)    //enforce rules with tagsj and tagsk
            {
	      if (enforce_rules[i].tagsk[k2]==k){
	        found = true;
	        break;
	      }
            }
            if (!found) (td->getA())[enforce_rules[i].tagi][j][k] = ZERO;
	  }	  
        }
      }
    }
  }
    
  //cerr<<"HMM2::applyrules enforce rules dones\n";
  // Normalize probabilities
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++){ 
      double sum=0;
      for(k=0; k<N; k++) 
        sum += (td->getA())[i][j][k];
      for(k=0; k<N; k++) {
        if (sum>0)
	  (td->getA())[i][j][k] = (td->getA())[i][j][k]/sum;
        else
	  (td->getA())[i][j][k] = 0;
      }
    }
  }
  //cerr<<"HMM2::applyrules Normalize probs done\n";

}

void 
HMM2::read_dictionary (FILE *fdic) {
  //cerr<<"HMM2::readDictionary\n";
  int i, k, nw=0;
  TaggerWord *word=NULL;
  set <TTag> tags;
  Collection &output = td->getOutput();
  
  MorphoStreamTrigram morpho_stream(fdic, true, td);
  
  // In the input dictionary there must be all punctuation marks, including the end-of-sentece mark
   
  word = morpho_stream.get_next_word();
  
  while (word) {
    if (++nw%10000==0) wcerr<<L'.'<<flush;
    
    tags = word->get_tags();

    if (tags.size()>0)
      k = output[tags];

    delete word;
    word = morpho_stream.get_next_word();
  }
  wcerr<<L"\n";
  
  // OPEN AMBIGUITY CLASS
  // It contains all tags that are not closed.
  // Unknown words are assigned the open ambiguity class
  k=output[td->getOpenClass()];

  int N = (td->getTagIndex()).size();  
  
  // Create ambiguity class holding one single tag for each tag.
  // If not created yet
  for(i = 0; i != N; i++) {
    set<TTag> amb_class;
    amb_class.insert(i);
    k=output[amb_class];
  }

  int M = output.size();
  
  wcerr<< N <<L" states and "<< M <<L" ambiguity classes\n";
  td->setProbabilities(N, M);
}

void
HMM2::filter_ambiguity_classes(FILE *in, FILE *out) {
  set<set<TTag> > ambiguity_classes;
  MorphoStreamTrigram morpho_stream(in, true, td);
  
  TaggerWord *word = morpho_stream.get_next_word();
  
  while(word) {
    set<TTag> tags = word->get_tags();

    if(tags.size() > 0) {     
      if(ambiguity_classes.find(tags) == ambiguity_classes.end()) {
	ambiguity_classes.insert(tags);
	word->outputOriginal(out);
	//wcerr<<word->get_string_tags()<<L"\n";
      }
    }
    delete word;
    word = morpho_stream.get_next_word();
  }
}

void 
HMM2::train (FILE *ftxt, int corpus_length, string savecountsfile) {
  //cerr<<"HMM2::BW\n";
  int i, j, k, k2, t, len, nw = 0;
  TaggerWord *word=NULL;
  TTag tag, pretag; 
  set<TTag> tags, pretags, prepretags;
  set<TTag>::iterator itag, jtag, ktag;
  map <int, double> gamma, gamma22;
  map <int, double>::iterator jt, kt;
  map <int, map <int, double> > phi;
  map <int, map <int, double> > gamma2;
  map < int, map <int,  map <int, double> > > alpha, beta, xsi2, phi2;
  map < int, map <int,  map <int, double> > >::iterator it;
  double prob, loli;              
  vector < set<TTag> > pending;
  Collection &output = td->getOutput();
  
  int ndesconocidas=0;
  map<int, double> ambclass_count;
  // alpha => forward probabilities
  // beta  => backward probabilities
  
  MorphoStreamTrigram morpho_stream(ftxt, true, td);

  loli = 0;
  tag = eos;
  tags.clear();
  tags.insert(tag);
  pending.push_back(tags);
  pending.push_back(tags);

  alpha[0].clear();      
  alpha[0][tag][tag] = 1;      //PROBLEM AREA: CHECK
  alpha[1][tag][tag] = 1;      //PROBLEM AREA: CHECK

  word = morpho_stream.get_next_word();

  while (word) {   
    //wcerr<<L"BW word: "<<word->get_superficial_form()<<L"\n";

    if (++nw%10000==0) wcerr<<L'.'<<flush;

    pretags = pending.back();
    prepretags = *(pending.end() - 2); // pending[pending.size()-2]

    tags = word->get_tags();    
    
    if (tags.size()==0) { // This is an unknown word
      tags = td->getOpenClass();
      ndesconocidas++;
    }
    
    if (output.has_not(tags)) {
      wstring errors;
      errors = L"A new ambiguity class was found. I cannot continue.\n";
      errors+= L"Word '"+word->get_superficial_form()+L"' not found in the dictionary.\n";
      errors+= L"New ambiguity class: "+word->get_string_tags()+L"\n";
      errors+= L"Take a look at the dictionary, then retrain.";
      fatal_error(errors);      
    }
    
    k = output[tags];    
    ambclass_count[k]++;
    len = pending.size();
    alpha[len].clear();     
      
    //Forward probabilities
    //wcerr<<L"BW:: forward prob start\n";
    for (itag=tags.begin(); itag!=tags.end(); itag++) {
      i=*itag;
      for (jtag=pretags.begin(); jtag!=pretags.end(); jtag++) {
         j=*jtag;
        for (ktag=prepretags.begin(); ktag!=prepretags.end(); ktag++) {
           k2=*ktag;
	   alpha[len][j][i] += alpha[len-1][k2][j]*(td->getA())[k2][j][i]*(td->getB())[j][i][k];
	   cerr<<alpha[len][j][i]<<"="<<alpha[len-1][k2][j]<<"x"<<(td->getA())[k2][j][i]<<"x"<<(td->getB())[j][i][k]<<"\n";
	}
      if (alpha[len][j][i]==0)
	{
          //alpha[len][j][i]=DBL_MIN;
	  cerr<<"DBL_MIN="<<alpha[len][j][i]<<"\n";
	}
      }
    }

    //wcerr<<L"BW:: forward prob end\n";
    if (!(tags.size()==1 && pretags.size()==1)) {
      pending.push_back(tags);
    } else {  // word and previous word are unambiguous
      //wcerr<<L"BW:: two unambiguous words\n";
      tag = *tags.begin(); 
      pretag = *pretags.begin(); 
      beta[0].clear();
      beta[0][pretag][tag] = 1;   
      prob = alpha[len][pretag][tag];         
      loli -= log(prob);  
      //wcerr<<L"prob="<<prob<<L" loli="<<loli<<L"\n";
      
      //wcerr<<L"BW:: backward prob start\n";
      for (t=0; t<len-1; t++) {  // loop from T-1 to 0	
        //wcerr<<L"t = "<<t<<L"\n";
        pretags = pending.back();
        prepretags = pending[pending.size()-2]; //*(pending.end() - 2)
	pending.pop_back();
   	k = output[tags];
	beta[1-t%2].clear();
	for (itag=tags.begin(); itag!=tags.end(); itag++) {
	  i=*itag;
          //wcerr<<"itag\n";
	  for (jtag=pretags.begin(); jtag!=pretags.end(); jtag++) {
	    j = *jtag;	      
            //wcerr<<"jtag\n";
	    for (ktag=prepretags.begin(); ktag!=prepretags.end(); ktag++) {
	      k2 = *ktag;	    
              //wcerr<<"A "<<k2<<" "<<j<<" "<<i<<" :"<<(td->getA())[k2][j][i]<<"\n";
              //wcerr<<"B "<<j<<" "<<i<<" "<<k<<" :"<<(td->getB())[j][i][k]<<"\n"; 
	      beta[1-t%2][k2][j] += (td->getA())[k2][j][i]*(td->getB())[j][i][k]*beta[t%2][j][i];
             // wcerr<<"... beta\n";
	      xsi2[k2][j][i] += alpha[len-t-1][k2][j]*(td->getA())[k2][j][i]*(td->getB())[j][i][k]*beta[t%2][j][i]/prob;
             // wcerr<<"... xsi2\n";
	    }
  	    gamma2[j][i] +=  alpha[len-t][j][i]*beta[t%2][j][i]/prob;	       
  	    gamma22[j] +=  alpha[len-t][j][i]*beta[t%2][j][i]/prob;	//sum of gamma2[j][i] over all i = tag_count 
  	    gamma[i] +=  alpha[len-t][j][i]*beta[t%2][j][i]/prob;	        //sum of phi[i][k] over all k = tag_count_for_emis
	    phi2[j][i][k] += alpha[len-t][j][i]*beta[t%2][j][i]/prob;               
	    phi[i][k] += alpha[len-t][j][i]*beta[t%2][j][i]/prob;               
            //wcerr<<"... gamma or phi\n";
	  }
	}
        tags=pretags;
      }
	
      //wcerr<<L"BW:: backward prob end\n";
      pretags.clear();
      pretags.insert(pretag);
      pending.push_back(pretags);
      tags.clear();
      tags.insert(tag);
      pending.push_back(tags);
      alpha[0].clear();
      alpha[0][pretag][tag] = 1;
    
      if((corpus_length>0)&&(nw>=corpus_length)&&(tag==eos)) {
        cerr<<nw<<" ";
        //End of training
        break;
      }
    }
    
    delete word; 
    word = morpho_stream.get_next_word();
  }  

  //wcerr<<L"BW:: out of loop\n";
  if ((pending.size()>1) || ((tag!=eos)&&(tag != (td->getTagIndex())[L"TAG_kEOF"])&&(pretag!=eos)&&(pretag != (td->getTagIndex())[L"TAG_kEOF"])) ) 
    wcerr<<L"Warning: Thee las tag is not the end-of-sentence-tag\n";
  
  
  if (savecountsfile!="") {
    cerr<<"Saving counts to file '"<<savecountsfile<<"' ... "<<flush;
    SmoothUtilsTrigram::save_counts(*td, savecountsfile, gamma22, gamma2, xsi2, ambclass_count, phi2, phi, gamma2, gamma);
    cerr<<"done.\n";
  }

  SmoothUtilsTrigram::calculate_smoothed_parameters(*td, gamma22, gamma2, xsi2, ambclass_count, phi2, phi, gamma2, gamma,  nw);

  wcerr<<L"Log="<<loli<<L"\n";
}

void 
HMM2::tagger(FILE *in, FILE *out, bool show_all_good_first) {
  //cerr<<"HMM2::Tagger\n";
  int i, j, k, k2, nw;
  TaggerWord *word=NULL;
  TTag tag, pretag, prepretag;
  
  set <TTag> tags, pretags, prepretags;
  set <TTag>::iterator itag, jtag, ktag;
  
  double prob, loli, x;
  int N = td->getN();  
  double alpha[2][N][N];
  vector<TTag> best[2][N][N];
  
  vector <TaggerWord> wpend; 
  int nwpend;
  
  MorphoStreamTrigram morpho_stream(in, debug, td);                             

  Collection &output = td->getOutput();
  
  loli = nw = 0;
   
  //Initialization
  tags.clear();
  pretags.clear();
  prepretags.clear();
  tags.insert(eos);
  pretags.insert(eos);
  alpha[0][eos][eos] = 1;       
   
  word = morpho_stream.get_next_word();
 
  while (word) {
   // wcerr<<L"HMM2::Tagger while loop start: word="<<word->get_superficial_form()<<L"\n";
    wpend.push_back(*word);    	    
    nwpend = wpend.size();
    
    prepretags = pretags; // Tags from the pre-previous word
    pretags = tags; // Tags from the previous word

    tags = word->get_tags();
  
    if (tags.size()==0) // This is an unknown word
      tags = td->getOpenClass();
                       
    if (output.has_not(tags)) {  // Encontrada una clase de ambigüedad desconocida hasta el momento      
      if (debug) {
        wstring errors;
	errors = L"A new ambiguity class was found. \n";
	errors+= L"Retraining the tagger is necessary so as to take it into account.\n";
	errors+= L"Word '"+word->get_superficial_form()+L"'.\n";
	errors+= L"New ambiguity class: "+word->get_string_tags()+L"\n";
	wcerr<<L"Error: "<<errors;
      }
      tags = find_similar_ambiguity_class(tags);
    } 
         
    k = output[tags];  //Ambiguity class the word belongs to
    
    //clear_array_double(alpha[nwpend%2], N);  
    //clear_array_vector(best[nwpend%2], N);
 
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	alpha[nwpend%2][i][j]=0.0;
        best[nwpend%2][i][j].clear();
      }
    }
    
    //cerr<<"HMM2::Tagger tags="<<tags.size()<<" pretags="<<pretags.size()<<" prepretags="<<prepretags.size()<<"\n";
    //Induction
    for (itag=tags.begin(); itag!=tags.end(); itag++) { //For all tag from the current word
      i=*itag;
      for (jtag=pretags.begin(); jtag!=pretags.end(); jtag++) {	//For all tags from the previous word
	j=*jtag;
        for (ktag=prepretags.begin(); ktag!=prepretags.end(); ktag++) {
           k2=*ktag;
	  x = alpha[1-nwpend%2][k2][j]*(td->getA())[k2][j][i]*(td->getB())[j][i][k];
          //cerr<<"HMM2::Tagger Induction x="<<x<<"\n";
          //cerr<<"HMM2::Tagger Induction alpha nwpend2 j i="<<alpha[nwpend%2][j][i]<<"\n";
	  //x = alpha[1-nwpend%2][j]*(td->getA())[j][i]*(td->getB())[i][k];
	  if (alpha[nwpend%2][j][i]<=x) {
          //cerr<<"HMM2::Tagger inside nwpend block, nwpend="<<nwpend<<"\n";
	    if (nwpend>1) 
	      best[nwpend%2][j][i] = best[1-nwpend%2][k2][j];
	    best[nwpend%2][j][i].push_back(i);
	    alpha[nwpend%2][j][i] = x;
	  }
        }
      }
    }
    //cerr<<"HMM2::Tagger alpha best done\n";
    
    //Backtracking
    if (tags.size()==1 && pretags.size()==1) {       
      cerr<<"HMM2::Tagger UNAMBIGUOUS pair!\n";
      tag = *tags.begin();      
      pretag = *pretags.begin();      
      
      prob = alpha[nwpend%2][pretag][tag];
      
      if (prob>0) 
	loli -= log(prob);
      else {
        if (debug)
	  wcerr<<L"Problem with word '"<<word->get_superficial_form()<<L"' "<<word->get_string_tags()<<L"\n";
      }
      //cerr<<"HMM2::Tagger print best pretag tag size "<<best[nwpend%2][pretag][tag].size()<<"\n";
      for (unsigned t=0; t<best[nwpend%2][pretag][tag].size(); t++) {
	if (show_all_good_first) {
          cerr<<"HMM2::Tagger print in backtracking good first\n";
	  wstring const &micad = wpend[t].get_all_chosen_tag_first(best[nwpend%2][pretag][tag][t], (td->getTagIndex())[L"TAG_kEOF"]);
	  fputws_unlocked(micad.c_str(), out); 
	} else { 
          cerr<<"HMM2::Tagger print in backtracking\n";
	  wpend[t].set_show_sf(show_sf);
	  wstring const &micad = wpend[t].get_lexical_form(best[nwpend%2][pretag][tag][t], (td->getTagIndex())[L"TAG_kEOF"]);
	  fputws_unlocked(micad.c_str(), out); 
	}
      }
      
      //Return to the initial state
      wpend.clear();   
      alpha[0][pretag][tag] = 1;
    }
    
    delete word;
    word = morpho_stream.get_next_word();    
  }
  
  if ((tags.size()>1)&&(debug)) {
    wstring errors;
    errors = L"The text to disambiguate has finished, but there are ambiguous words that has not been disambiguated.\n";
    errors+= L"This message should never appears. If you are reading this ..... these are very bad news.\n";
    wcerr<<L"Error: "<<errors;
  }  
}


void
HMM2::print_A() {
  int i,j,k;
    
  cout<<"TRANSITION MATRIX (A)\n------------------------------\n";  
  for(i=0; i != td->getN(); i++)
    for(j=0; j != td->getN(); j++) {
      for(k=0; k != td->getN(); k++) {
        cout<<"A["<<i<<"]["<<j<<"]["<<k<<"] = "<<(td->getA())[i][j][k]<<"\n";
      }    
    }    
}

void
HMM2::print_B() {
  int i,j,k;  

  cout<<"EMISSION MATRIX (B)\n-------------------------------\n";
  for(i=0; i != td->getN(); i++)
    for(j=0; j != td->getN(); j++){
      for(k=0; k != td->getM(); k++) {
        Collection &output = td->getOutput();
        if(output[k].find(j)!=output[k].end())
          cout<<"B["<<i<<"]["<<j<<"]["<<k<<"] = "<<(td->getB())[i][j][k]<<"\n";
      }
    }
}

void 
HMM2::print_ambiguity_classes() {
  set<TTag> ambiguity_class;
  set<TTag>::iterator itag;
  cout<<"AMBIGUITY CLASSES\n-------------------------------\n";
  for(int i=0; i != td->getM(); i++) {
    ambiguity_class = (td->getOutput())[i];
    cout <<i<<": ";
    for (itag=ambiguity_class.begin(); itag!=ambiguity_class.end(); itag++) {
      cout << *itag <<" ";
    }
    cout << "\n";
  }
}   

set<TTag>
HMM2::find_similar_ambiguity_class(set<TTag> c) {
  int size_ret = -1;
  set<TTag> ret=td->getOpenClass(); //Se devolverá si no encontramos ninguna clase mejor
  bool skeep_class;
  Collection &output = td->getOutput();

  for(int k=0; k<td->getM(); k++) {
    if ((((int)output[k].size())>((int)size_ret)) && (((int)output[k].size())<((int)c.size()))) {
      skeep_class=false;
      // Test if output[k] is a subset of class
      for(set<TTag>::iterator it=output[k].begin(); it!=output[k].end(); it++) {
        if (c.find(*it)==c.end()) { 
	   skeep_class=true; //output[k] is not a subset of class
	   break;
	}
      }
      if (!skeep_class) {
        size_ret = output[k].size();
	     ret = output[k];
      }
    }
  }
  return ret;
}

double 
HMM2::log_add (double x, double y) {
  if (x < y) {
    double temp;
    temp = x;
    x = y;
    y = temp;
  }
  if (y == -numeric_limits<double>::infinity()) {  //log(0) = - Infinity
    cerr<<"log_add: one of the log probs is null"<<endl;
    return x;
  } else {
    double diff = y - x;
    return x + log(1.0 + exp(diff));
  }
}
