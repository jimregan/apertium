#include <apertium/perceptron_tagger.h>
#include <apertium/perceptron_spec.h>
#include <apertium/wchar_t_exception.h>
#include <algorithm>
#include <map>
#include <set>

namespace Apertium {

/*
double PerceptronTagger::local_weight(std::set<std::string> &feats) {
  std::set<std::string>::const_iterator feat_iter = feats.begin();
  double out = 0;
  for (;feat_iter!=feats.end();feat_iter++) {
    out += tdpercep.weights[*feat_iter];
  }
  return out;
}
*/

PerceptronTagger::PerceptronTagger(basic_Tagger::Flags flags) : basic_Tagger(flags) {};

PerceptronTagger::~PerceptronTagger() {};

void PerceptronTagger::tag(Stream &in, std::wostream &out) const {
  SentenceStream::SentenceTagger::tag(in, out);
}

TaggedSentence
PerceptronTagger::tagSentence(const Sentence &untagged_sent) const {
  const size_t sent_len = untagged_sent.size();

  std::vector<AgendaItem> agenda;
  agenda.reserve(spec.beam_width);
  agenda.push_back(AgendaItem());
  agenda.back().tagged.reserve(sent_len);

  std::vector<std::string> feat_vec_delta;
  std::vector<Analysis>::const_iterator analys_it;
  std::vector<AgendaItem>::const_iterator agenda_it;
  std::vector<Morpheme>::const_iterator wordoid_it;

  for (size_t token_idx = 0; token_idx < sent_len; token_idx++) {
    const std::vector<Analysis> &analyses =
        untagged_sent[token_idx].TheLexicalUnit->TheAnalyses;

    std::vector<AgendaItem> new_agenda;
    new_agenda.reserve(spec.beam_width * analyses.size());

    if (analyses.size() == 1) {
      extendAgendaAll(agenda, analyses.front());
      continue;
    } else {
      extendAgendaAll(agenda, Optional<Analysis>());
      continue;
    }

    for (agenda_it = agenda.begin(); agenda_it != agenda.end(); agenda_it++) {
      for (analys_it = analyses.begin(); analys_it != analyses.end(); analys_it++) {
        const std::vector<Morpheme> &wordoids = analys_it->TheMorphemes;

        new_agenda.push_back(*agenda_it);
        AgendaItem &new_agenda_item = new_agenda.back();
        new_agenda_item.tagged.push_back(*analys_it);

        for (wordoid_it = wordoids.begin(); wordoid_it != wordoids.end(); wordoid_it++) {
          int wordoid_idx = wordoid_it - wordoids.begin();
          feat_vec_delta.clear();
          spec.get_features(new_agenda_item.tagged, untagged_sent,
                            token_idx, wordoid_idx, feat_vec_delta);
          new_agenda_item.score += tdpercep.weights * feat_vec_delta;
        }
      }
    }
    // Apply the beam
    size_t new_agenda_size = std::min((size_t)spec.beam_width, new_agenda.size());
    agenda.resize(new_agenda_size);
    std::partial_sort_copy(new_agenda.begin(), new_agenda.end(),
                           agenda.begin(), agenda.end());
  }

  return agenda.front().tagged;
}

void PerceptronTagger::outputLexicalUnit(
    const LexicalUnit &lexical_unit, const Optional<Analysis> analysis,
    std::wostream &output) const {
  StreamTagger::outputLexicalUnit(lexical_unit, analysis, output);
}

void PerceptronTagger::trainSentence(const TrainingSentence &sentence)
{
  const TaggedSentence &tagged_sent = sentence.first;
  const Sentence &untagged_sent = sentence.second;
  assert(tagged_sent.size() == untagged_sent.size());
  const size_t sent_len = tagged_sent.size();

  std::vector<TrainingAgendaItem> agenda;
  agenda.reserve(spec.beam_width);
  agenda.push_back(TrainingAgendaItem());
  agenda.back().tagged.reserve(sent_len);
  std::vector<TrainingAgendaItem>::const_iterator correct_agenda_it
      = agenda.begin();

  TrainingAgendaItem correct_sentence;
  correct_sentence.tagged.reserve(sent_len);

  std::vector<std::string> feat_vec_delta;
  std::vector<Analysis>::const_iterator analys_it;
  std::vector<TrainingAgendaItem>::const_iterator agenda_it;
  std::vector<Morpheme>::const_iterator wordoid_it;

  for (size_t token_idx = 0; token_idx < sent_len; token_idx++) {
    const TaggedToken &tagged_tok(tagged_sent[token_idx]);
    const StreamedType &untagged_tok(untagged_sent[token_idx]);
    correct_sentence.tagged.push_back(tagged_tok);

    const std::vector<Analysis> &analyses =
        untagged_tok.TheLexicalUnit->TheAnalyses;

    std::vector<TrainingAgendaItem> new_agenda;
    new_agenda.reserve(spec.beam_width * analyses.size());

    if (analyses.size() <= 1 || !tagged_tok) {
      // Case |analyses| = 0, nothing we can do
      // Case !tagged_tok, |analyses| > 0, no point penalising a guess which
      //   can only be incorrect when there's no correct answer
      // Case |analyses| = 1, everything will cancel out anyway
      if (analyses.size() == 1) {
        extendAgendaAll(agenda, analyses.front());
        continue;
      } else {
        extendAgendaAll(agenda, Optional<Analysis>());
        continue;
      }
    }

    bool correct_available = false;
    for (agenda_it = agenda.begin(); agenda_it != agenda.end(); agenda_it++) {
      for (analys_it = analyses.begin(); analys_it != analyses.end(); analys_it++) {
        const std::vector<Morpheme> &wordoids = analys_it->TheMorphemes;

        new_agenda.push_back(*agenda_it);
        TrainingAgendaItem &new_agenda_item = new_agenda.back();
        new_agenda_item.tagged.push_back(*analys_it);

        for (wordoid_it = wordoids.begin(); wordoid_it != wordoids.end(); wordoid_it++) {
          int wordoid_idx = wordoid_it - wordoids.begin();
          feat_vec_delta.clear();
          spec.get_features(new_agenda_item.tagged, untagged_sent,
                            token_idx, wordoid_idx, feat_vec_delta);
          new_agenda_item.vec += feat_vec_delta;
          new_agenda_item.score += tdpercep.weights * feat_vec_delta;
          if (agenda_it == correct_agenda_it && *analys_it == *tagged_tok) {
            correct_sentence = new_agenda_item;
            correct_available = true;
          }
        }
      }
    }
    if (!correct_available) {
      // XXX: Should possibly allow a flag to convert this error to a warning
      // and simply skip the sentence at this point
      std::wstringstream what_;
      what_ << L"Tagged analysis unavailable in untagged/ambigous input.\n";
      what_ << L"Available:\n";
      for (analys_it = analyses.begin(); analys_it != analyses.end(); analys_it++) {
        what_ << *analys_it << L"\n";
      }
      what_ << L"Required: " << *tagged_tok << L"\n";
      throw Apertium::wchar_t_Exception::PerceptronTagger::CorrectAnalysisUnavailable(what_);
    }
    // Apply the beam
    size_t new_agenda_size = std::min((size_t)spec.beam_width, new_agenda.size());
    agenda.resize(new_agenda_size);
    std::partial_sort_copy(new_agenda.begin(), new_agenda.end(),
                           agenda.begin(), agenda.end());

    // Early update "fallen off the beam"
    bool any_match = false;
    for (agenda_it = agenda.begin(); agenda_it != agenda.end(); agenda_it++) {
      if (agenda_it->tagged == correct_sentence.tagged) {
        correct_agenda_it = agenda_it;
        any_match = true;
        break;
      }
    }
    if (!any_match) {
      tdpercep.weights -= agenda.front().vec;
      tdpercep.weights += correct_sentence.vec;
      return;
    }
  }
  // Normal update
  if (agenda.front().tagged != correct_sentence.tagged) {
    tdpercep.weights -= agenda.front().vec;
    tdpercep.weights += correct_sentence.vec;
  }
}

void PerceptronTagger::train(Stream&) {}  // dummy

void PerceptronTagger::train(
    Stream &tagged,
    Stream &untagged,
    int iterations) {
  std::map<std::string, unsigned int> tstamps;
  FeatureVec totals;

  TrainingCorpus tc(tagged, untagged);
  for (int i=0;i<iterations;i++) {
    tc.shuffle();
    std::vector<TrainingSentence>::const_iterator si;
    for (si = tc.sentences.begin(); si != tc.sentences.end(); si++) {
      trainSentence(*si);
    }
  }
}

void PerceptronTagger::serialise(std::ostream &serialised) const {};

void PerceptronTagger::deserialise(std::istream &serialised) {} // dummy

void PerceptronTagger::deserialise(TaggerData const&) {} // dummy

template <typename T> void
PerceptronTagger::extendAgendaAll(
    std::vector<T> &agenda,
    Optional<Analysis> analy) {
  typename std::vector<T>::iterator agenda_it;
  for (agenda_it = agenda.begin(); agenda_it != agenda.end(); agenda_it++) {
    agenda_it->tagged.push_back(analy);
  }
}

bool operator<(PerceptronTagger::AgendaItem &a,
               PerceptronTagger::AgendaItem &b) {
  return a.score < b.score;
};
}
