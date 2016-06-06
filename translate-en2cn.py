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


# word is an annotation
def translate_one_word_deprecated(word):
   word = word.strip()
   debug_untranslated = ''
   if word == '':
      return ''
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
      num_translated = 0
      lst = word.split()
      n = len(lst)      
      for i in range(n): # number of passes in translating
         translate_len = []
         translate_lst = []
         lst = word.split()
         n = len(lst)
         found = 0
         for j in range(n+1,i,-1):
            s = ' '.join(lst[i:j])
            if s in d:
               found = 1
               translate_len.append(j-i)
               translate_lst.append(' '.join(lst[:i]) + ' ' + d[s] + ' ' + ' '.join(lst[j:]))
               break
         max_len = 0
         if len(translate_lst) > 0:
            max_len = max(translate_len)
            index = translate_len.index(max_len)
            word =  translate_lst[index]
            num_translated += max_len
         else:
            debug_untranslated += ' ' + ' '.join(lst[i:j])
         if num_translated == n:
            break
      #return debug_untranslated.strip()
      return word.strip()


# word is an annotation
def translate_one_word(word):
   word = word.strip()
   if word == '':
      return ''
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
            translate += ' ' + lst[i]
            i += 1
         else:
            translate += ' ' + d[s]
            i += j - i;
      return translate.strip()


# word is an annotation
def find_untranslable_words(word):
   word = word.strip()
   if word == '':
      return ''
   i = word.find('gene:')
   if  i > -1:
      return ''
   i = word.find('locus:')
   if i > -1:
      return ''
   i = word.find('Publication:')
   if i > -1:
      return ''
   i = word.find('AnalysisReference:')
   if i > -1:
      return ''
   i = word.find('Communication:')
   if i > -1:
      return '' 
   
   if word in d:
      return ''
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
            translate += ' ' + lst[i]
            i += 1
         else:
            i += j - i;
      result = translate.strip()
      if result.startswith('GO:') or len(result) == 3 or result.startswith('AT') or result.isnumeric() \
         or result.startswith('20') or result.startswith('TAIR') or result.startswith('TIGR') \
         or result.find(':') > 0:
         return ''
      return result


def unique_words_to_files(words, file):
   
   lst = words.strip().split()
   d = {}
   for x in lst:
      if x in d:
         d[x] += 1
      else:
         d[x] = 0
   sorted_d = sorted(d.items(), key=operator.itemgetter(1), reverse=True)
   for (k, v) in sorted_d:
      file.write('     \'' + k + '\'' + ':' + '\'' + str(v) + '\',' + '\n')
      
# Write translation to file
tfile = codecs.open('arabidopsis-GOSLIM-translation.txt','w','utf-8')
untranslated_str = ''
content = file_to_lst('ATH_GO_GOSLIM.txt')
gene_name_d = make_gene_name_dict('gene_aliases_20140331.txt')
for line in content:
   lst = line.strip().split('\t')
   gene = lst[0]
   #gene_name = gene_name_d.get(gene, '?')
   description = gene + '\t'
   for s in lst[1:]:
      t = translate_one_word(s);
      description += '\t' + t.replace(' ', '')
      #u = find_untranslable_words(s)
      #untranslated_str += ' ' + u
   tfile.write(description + '\n')
tfile.close()

