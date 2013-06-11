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
 *  Sliding-Window Part of Speech Tagger (SWPoST) implementation (source)
 *
 *  @author   Gang Chen - pkuchengang@gmail.com
 */


#include <apertium/swpost.h>
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

SWPoST::SWPoST(TaggerData *t)
{
  this->td = t;
}

SWPoST::~SWPoST()
{
}

void
SWPoST::init()
{
}

void SWPoST::init_probabilities(FILE *ftxt) {
	int N = td->getN();
	int M = td->getM();
	int i, j, k, s_left, s_right, nw = 0; //sigma_left, sigma_right

#ifdef __GNUC__
	float para_matrix[M][M][N]; //M = Number of ambiguity classes, N = Number of fine tags
	int para_matrix_sum[M][M];
	for (i = 0; i < M; ++i) {
		for (j = 0; j < M; ++j) {
			para_matrix_sum[i][j] = 0;
			for (k = 0; k < N; ++k ) {
				para_matrix[i][j][k] = 0;
			}
		}
	}
#else
	vector<vector<vector< float > > > para_matrix(M, vector<vector<float> >(M, vector<float>(N, 0)));
	vector<vector< int > > para_matrix_sum(M, vector<int>(M, 0));
#endif

	Collection &output = td->getOutput();

	MorphoStream lexmorfo(ftxt, true, td);
	
	set<TTag> tags_left, tags, tags_right;

	TaggerWord *word_left = lexmorfo.get_next_word();
	TaggerWord *word = lexmorfo.get_next_word();
	TaggerWord *word_right = lexmorfo.get_next_word();
	if (word_left == NULL || word == NULL || word_right == NULL) {
		return;
	}

	//We count each element of the para matrix
	while (word_right != NULL) {

		if (++nw % 10000 == 0)
			wcerr << L'.' << flush;

		tags_left = word_left->get_tags();
		if (tags_left.size() == 0) { //This is an unknown word
			tags_left = td->getOpenClass();
		}
		tags = word->get_tags();
		if (tags.size() == 0) { //This is an unknown word
			tags = td->getOpenClass();
		}
		tags_right = word_right->get_tags();
		if (tags_right.size() == 0) { //This is an unknown word
			tags_right = td->getOpenClass();
		}

		if (output.has_not(tags) || output.has_not(tags_left)
				|| output.has_not(tags_right)) {
			wstring errors;
			errors = L"A new ambiguity class was found. I cannot continue.\n";
			errors += L"Word '" + word->get_superficial_form() + L"' not found in the dictionary.\n";
			errors += L"New ambiguity class: " + word->get_string_tags() + L"\n";
			errors += L"Take a look at the dictionary and at the training corpus. Then, retrain.";
			fatal_error(errors);
		}

		s_left = output[tags_left];
		s_right = output[tags_right];

		for (set<TTag>::iterator iter = tags.begin();
				iter != tags.end(); ++iter) {
			para_matrix[s_left][s_right][*iter]++;
		}
		para_matrix_sum[s_left][s_right]++;

		delete word_left;
		word_left = word;
		word = word_right;
		word_right = lexmorfo.get_next_word();
	}

	// get normalized count, by context of "s_left ... s_right"
	for (i = 0; i < M; ++i) {
		for (j = 0; j < M; ++j) {
			for (k = 0; k < N; ++k) {
				para_matrix[i][j][k] = para_matrix[i][j][k]
						/ para_matrix_sum[i][j];
			}
		}
	}
	//td->setSWPoSTProbabilities(N, M, para_matrix);
	for (i = 0; i < M; ++i) {
        	for (j = 0; j < M; ++j) {
                	for (k = 0; k < N; ++k) {
				td->getC()[i][j][k] = para_matrix[i][j][k];
			}
		}
	}

	wcerr << L"\n";
}

void SWPoST::read_dictionary(FILE *fdic) {
	int i, k, nw = 0;
	TaggerWord *word = NULL;
	set < TTag > tags;
	Collection &output = td->getOutput();

	MorphoStream morpho_stream(fdic, true, td);

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
	k = output[td->getOpenClass()];

	int N = (td->getTagIndex()).size();

	// Create ambiguity class holding one single tag for each tag.
	// If not created yet
	for (i = 0; i != N; i++) {
		set < TTag > amb_class;
		amb_class.insert(i);
		k = output[amb_class];
	}

	int M = output.size();

	wcerr << N << L" states and " << M << L" ambiguity classes\n";
	td->setSWPoSTProbabilities(N, M);
}

void SWPoST::train(FILE *ftxt) {
	int N = td->getN();
	int M = td->getM();
	int i, j, k, s_left, s_right, nw = 0; //sigma_left, sigma_right

#ifdef __GNUC__
	float para_matrix_new[M][M][N]; //M = Number of ambiguity classes, N = Number of tags (states)
	for (i = 0; i < M; ++i) {
		for (j = 0; j < M; ++j) {
			for (k = 0; k < N; ++k ) {
				para_matrix_new[i][j][k] = 0;
			}
		}
	}
#else
	vector<vector<vector<float> > > para_matrix_new(M, vector<vector<float > >(M, vector<float>(N, 0)));
#endif

	set<TTag> tags_left, tags, tags_right;

	Collection &output = td->getOutput();

	MorphoStream morpho_stream(ftxt, true, td);

	TaggerWord *word_left = morpho_stream.get_next_word();
	TaggerWord *word = morpho_stream.get_next_word();
	TaggerWord *word_right = morpho_stream.get_next_word();
	if (word_left == NULL || word == NULL || word_right == NULL) {
		return;
	}

	while (word_right) {

		//wcerr<<L"Enter para continuar\n";
		//getchar();
		if (++nw % 10000 == 0)
			wcerr << L'.' << flush;

		tags_left = word_left->get_tags();
		if (tags_left.size() == 0) { //This is an unknown word
			tags_left = td->getOpenClass();
		}
		tags = word->get_tags();
		if (tags.size() == 0) { //This is an unknown word
			tags = td->getOpenClass();
		}
		tags_right = word_right->get_tags();
		if (tags_right.size() == 0) { //This is an unknown word
			tags_right = td->getOpenClass();
		}

		if (output.has_not(tags_left) || output.has_not(tags)
				|| output.has_not(tags_right)) {
			wstring errors;
			errors = L"A new ambiguity class was found. I cannot continue.\n";
			errors += L"Word '" + word->get_superficial_form() + L"' not found in the dictionary.\n";
			errors += L"New ambiguity class: " + word->get_string_tags() + L"\n";
			errors += L"Take a look at the dictionary and at the training corpus. Then, retrain.";
			fatal_error(errors);
		}

		s_left = output[tags_left];
		s_right = output[tags_right];

		double normalization = 0;
		for (set<TTag>::iterator iter = tags.begin();
				iter != tags.end(); ++iter) {
			normalization += (td->getC())[s_left][s_right][*iter];
		}
		for (set<TTag>::iterator iter = tags.begin();
				iter != tags.end(); ++iter) {
			para_matrix_new[s_left][s_right][*iter] += 1.0 / normalization;
		}

		delete word_left;
		word_left = word;
		word = word_right;
		word_right = morpho_stream.get_next_word();
	}

	//td-setSWPoSTProbabilities(N, M, (double ***)para_matrix_new);
	for (i = 0; i < M; ++i) {
        	for (j = 0; j < M; ++j) {
                	for (k = 0; k < N; ++k) {
				td->getC()[i][j][k] = para_matrix_new[i][j][k];
			}
		}
	}

	wcerr << L"\n";
}

void
SWPoST::print_para_matrix() {
}


//TODO
void 
SWPoST::taggerSWPoST(FILE *in, FILE *out, bool show_all_good_first) {
wcerr << L"xxx taggerSWPoST:" << endl;
  int s_left, s_right;
  TaggerWord *word_left = NULL;
  TaggerWord *word = NULL;
  TaggerWord *word_right = NULL;
  
  set <TTag> tags_left, tags, tags_right;
  
  int M = td->getM();  
  int N = td->getN();  
wcerr << L"M, N =" << M << L"," << N << endl;
  MorphoStream morpho_stream(in, true, td);                             
  
  Collection &output = td->getOutput();
  
  word_left = morpho_stream.get_next_word(); 
  word = morpho_stream.get_next_word();
  word_right = morpho_stream.get_next_word();
wcerr << L"xxx: " << *word_left << endl;
wcerr << L"xxx: " << *word << endl;
wcerr << L"xxx: " << *word_right << endl;
  if (word_left == NULL || word == NULL || word_right == NULL) {
      return;
  }

  wstring micad;

  //first word, naive.
  micad = word_left->get_lexical_form((TTag &)*(word_left->get_tags().begin()), (td->getTagIndex())[L"TAG_kEOF"]);
  fputws_unlocked(micad.c_str(), out);
 
  while (word_right) {

	wcerr << L"xxx: word_right=" << *word_right << endl;

  	tags_left = word_left->get_tags();
	if (tags_left.size() == 0) {
		tags_left = td->getOpenClass();
	}
	tags = word->get_tags();
	if (tags.size() == 0) {
		tags = td->getOpenClass();
	}
	tags_right = word_right->get_tags();
	if (tags_right.size() == 0) {
		tags_right = td->getOpenClass();
	}

	if (output.has_not(tags_left) || output.has_not(tags) || output.has_not(tags_right)) {
		wstring errors;
		errors = L"A new ambiguity class was found. \n";
		errors+= L"Retraining the tagger is necessary so as to take it into account.\n";
		errors+= L"Word '"+word->get_superficial_form()+L"'.\n";
		errors+= L"New ambiguity class: "+word->get_string_tags()+L"\n";
		wcerr<<L"Error: "<<errors;
	}
		
	s_left =  output[tags_left];
	s_right = output[tags_right];

	double max = 0;
	TTag tag_max = 0;
	for (set<TTag>::iterator iter = tags.begin(); iter != tags.end(); ++iter) {
		double n = (td->getC())[s_left][s_right][*iter];
		if (n > max) {
			max = n;
			tag_max = *iter;
		}
	}
	micad = word->get_lexical_form(tag_max, (td->getTagIndex())[L"TAG_kEOF"]);
	fputws_unlocked(micad.c_str(), out);

        word_left = word;
	word = word_right;
	word_right = morpho_stream.get_next_word();
  }

  //last word, naive.
  micad = word_right->get_lexical_form((TTag&)*(word_right->get_tags().begin()), (td->getTagIndex())[L"TAG_kEOF"]);
  fputws_unlocked(micad.c_str(), out);

  wcerr << L"\n";
}

