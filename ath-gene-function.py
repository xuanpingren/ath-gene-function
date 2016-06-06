#!/usr/bin/python
# -*- coding: utf-8 -*
import cgi
import cgitb
cgitb.enable()
import random
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


def convert_to_locus(gene, filename):
   gene = gene.strip().lower()
   if gene == '':
      return 'NONE'
   f = open(filename)
   for line in f:
      lst = line.split('\t')
      locus = lst[0]
      short_name = lst[1].strip().lower()
      short_name_lst = short_name.split()
      long_name = lst[2].strip().lower()
      long_name_lst = long_name.split()
      gene_lst = gene.split()
      n = len(gene_lst)
      count1 = 0
      count2 = 0
      for x in gene_lst:
         if x in short_name_lst:
            count1 += 1
         if x in long_name_lst:
            count2 += 1
      if count1 == n or count2 == n:
         return locus
         
   return 'NONE'


def get_gene_name(gene, file):

   name = ''
   f = open(file)
   count = 0
   for line in f:
      lst = line.split('\t')
      if gene == lst[0]:
         name += ' '.join(lst[1:]).strip() + ';'
         count += 1
      if count > 10:
         break
   f.close()
   return name


# word is an annotation
def translate_one_word(word):
   if word in d:
      return d[word]
   else:
      #for x in word.split():
      #   x = x.strip()
      #   if x in d:
      #      t += d[x] +  ' '
      #   else:
      #      t += x + ' '
      for k in range(5): # number of passes in translating

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

def find_gene_annotation(gene, lst):
   
   p = []
   c = []
   f = []
   for x in lst:
      L = x.strip().split('\t')
      t = L[7]
      a = L[4]
      if L[0] == gene:
         if t  == 'P' and not a in p:
            p.append(a)
         elif t == 'C' and not a in c:
            c.append(a)
         elif t == 'F' and not a in f:
            f.append(a)
   s = ''
   if p:
      s += '<br/><b>生物过程</b><br/>' + '<br/>'.join(translate_words(p)) + '<br/>'
   if f:
      s += '<br/><b>分子功能</b><br/>' + '<br/>'.join(translate_words(f)) + '<br/>'
   if c:
      s += '<br/><b>所在位置</b><br/>' + '<br/>'.join(translate_words(c)) + '<br/>'
   return s


def find_genes_with_description(desc, lst):

   result_d = {}
   desc_lst = desc.split()
   for x in lst:
      L = x.strip().split('\t')
      t = L[7]
      a = L[4]
      a_lst = a.split()
      a_lst = [x.lower() for x in a_lst] # convert each element to lower case
      g = L[0]
      count = 0
      for y in desc_lst:
         if y in a_lst:
            count += 1
      if count == len(desc_lst):
         t = result_d.get(g, '')
         if t.find(a) == -1:
            t +=  a + '; '
         result_d[g] = t 

   return result_d

def in_locus(gene, content):
   for x in content:
      if gene == x.split('\t')[0]:
         return True
   return False

print("Content-type: text/html; charset=utf-8")
print("")

form = cgi.FieldStorage()
gene = ''
description = ''
if 'gene_name' in form:
   gene = form['gene_name'].value
   gene = gene.strip().upper()
if 'gene_function' in form:
   description = form['gene_function'].value
   description = description.strip().lower()
content = file_to_lst('ATH_GO_GOSLIM.txt')
invalid_gene = False
if gene == '':
   invalid_gene = True
if invalid_gene  and description == '':
   print("无效基因 %s<br/>" % gene)
   exit()
if not in_locus(gene, content):
   locus = convert_to_locus(gene, 'gene_aliases_20140331.txt')
   if locus != 'NONE':
      gene = locus
   else:
      invalid_gene = True



result = ''
if not invalid_gene:
   result = find_gene_annotation(gene, content)
   gene_name = get_gene_name(gene, 'gene_aliases_20140331.txt')

print("<!DOCTYPE html>")
print("<html>")
print("<head><title>%s  %s</title></head>" % (gene, description))
print("<body>")
if result != '':
   print("<p>基因：<font color=\"blue\">%s</font><br/>名字：%s</p>" % (gene, gene_name))
   print("%s" % result)
   print("<br/>")    
elif not invalid_gene:
   print("找不到。")


if description != '':
   gene_annot_d = find_genes_with_description(description, content)
   if gene_annot_d:
      print("<p><b>有%s功能的基因</b></p>" % (description))
   else:
      exit()
   print("<ul>")
   d = make_gene_name_dict('gene_aliases_20140331.txt')
   for x in gene_annot_d:
      g = x
      a = gene_annot_d[g]
      gn = d.get(g, '?')
      print("<li><font color=\"blue\">%s</font> (%s) <br/> %s</li><br/>" % (g, gn, a))
   print("<ul>")
print("</body>")
print("</html>")
      
      
