#!/usr/bin/python
# -*- coding: utf-8 -*
import random
import codecs
from word_dictionary import *

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


# word is an annotation
def translate_one_word(word):
   i = word.find('gene:')
   if  i > -1:
      return '基因' + word[4:]
   i = word.find('locus:')
   if i > -1:
      return '基因座' + word[5:]
   i = word.find('Publication:')
   if i > -1:
      return '出版物' + word[11:]
   i = word.find('AnalysisReference:')
   if i > -1:
      return '分析参考' + word[17:]
   i = word.find('Communication:')
   if i > -1:
      return '交流' + word[13:]   
   
   if word in d:
      return d[word]
   else:
      for k in range(8): # number of passes in translating

         translate_len = []
         translate_lst = []
         lst = word.split()
         n = len(lst)
         for i in range(n):
            found = 0
            for j in range(n,i,-1):
               s = ' '.join(lst[i:j])
               if s in d:
                  found = 1
                  translate_len.append(j-i)
                  translate_lst.append(' '.join(lst[:i]) + ' ' + d[s] + ' ' + ' '.join(lst[j:]))
                  break
         if len(translate_lst) > 0:
            index = translate_len.index(max(translate_len))
            word =  translate_lst[index]
         else:
            return word

      return word

def translate_words(word_lst):
   lst = []
   for x in word_lst:
      lst.append(translate_one_word(x))
   return lst


# Write translation to file
tfile = codecs.open('arabidopsis-GOSLIM-translation.txt','w','utf-8')
content = file_to_lst('ATH_GO_GOSLIM.txt')
gene_name_d = make_gene_name_dict('gene_aliases_20140331.txt')
for line in content[:100]:
   lst = line.strip().split('\t')
   gene = lst[0]
   #gene_name = gene_name_d.get(gene, '?')
   description = gene + '\t'
   for s in lst[1:]:
      t = translate_one_word(s);
      description += '\t' + t.replace(' ', '')
   tfile.write(description + '\n')
tfile.close()
