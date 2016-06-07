#!/usr/bin/python
# -*- coding: utf-8 -*
import operator
import codecs
from word_dictionary import d

def file_to_lst(filename):
   f = open(filename)
   lst = f.readlines()
   f.close()
   return lst


def make_gene_name_dict(file):
   d = {}
   f = open(file)
   for line in f:
      lst = line.split('\t')
      gene = lst[0]
      name = d.get(gene, '')
      name += ' '.join(lst[1:]).strip() + ';'
      d[gene] = name
   f.close()
   return d


def split_and_translate_component(w, d):
   lst = w.split('-')
   result = ''
   for x in lst:
      if x in d:
         result += d[x] + '-'
      else:
         result += x + '-'
   return result[:-1]
   
# word is an annotation
def translate_one_word(word, d):
   word = word.strip()
   if word == '':
      return ''
   
   d0 = {'gene:':'基因', 'locus:':'基因座', 'Publication:':'出版物', 'AnalysisReference:':'分析参考', 'Communication:':'交流'}
   for k in d0:
      i = word.find(k)
      if i > -1:
         length = len(k)
         return d0[k] + word[length-1:]

   if word in d:
      return d[word]
   else:
      i = 0
      translate = ''
      lst = word.split()
      n = len(lst)      
      while i < n: # number of passes in translating
         found = 0
         for j in range(n+1,i,-1):
            s = ' '.join(lst[i:j])
            if s in d:
               found = 1
               break
         if found == 0:
            translate += ' ' + split_and_translate_component(lst[i], d)
            i += 1
         else:
            translate += ' ' + d[s]
            i += j - i;
      return translate.strip()



# Write translation to file
tfile = codecs.open('arabidopsis-GOSLIM-translation.txt','w','utf-8')
untranslated_str = ''
content = file_to_lst('ATH_GO_GOSLIM.txt')
gene_name_d = make_gene_name_dict('gene_aliases_20140331.txt')
for line in content:
   lst = line.strip().split('\t')
   gene = lst[0]
   description = gene + '\t'
   for s in lst[1:]:
      t = translate_one_word(s, d);
      description += '\t' + t.replace(' ', '')
   tfile.write(description + '\n')
tfile.close()

