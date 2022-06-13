rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)
library(sampleSelection) # Atenção!

# Abre a base de dados "mroz" do Wooldridge.
# T. A. Mroz (1987), “The Sensitivity of an Empirical Model of Married Women’s Hours of Work to Economic and Statistical Assumptions,” Econometrica 55, 765-799.
data(mroz)

# inlf: =1 if in lab frce, 1975
# hours: hours worked, 1975
# kidslt6: # kids < 6 years
# kidsge6: # kids 6-18
# age: woman’s age in yrs
# educ: years of schooling
# wage: est. wage from earn, hrs
# repwage: rep. wage at interview in 1976
# hushrs: hours worked by husband, 1975
# husage: husband’s age
# huseduc: husband’s years of schooling
# huswage: husband’s hourly wage, 1975
# faminc: family income, 1975
# mtr: fed. marg. tax rte facing woman
# motheduc: mother’s years of schooling
# fatheduc: father’s years of schooling
# unem: unem. rate in county of resid.
# city: =1 if live in SMSA
# exper: actual labor mkt exper
# nwifeinc: (faminc - wage*hours)/1000
# lwage: log(wage)
# expersq: exper^2

### SUMÁRIO ###
# 1 -- ESTIMAÇÃO (TRUNCAGEM INDICENTAL)

# ----------------------------------------------------------------------
# MODELO DE SELEÇÃO -- ESTIMAÇÃO (TRUNCAGEM INDICENTAL)
# ----------------------------------------------------------------------
# Nos casos em que é desejável aplicar o "método Heckit", não temos dificuldades adicionais de interpretação como nos outros modelos do capítulo. A implementação no R é automática e não é recomendável fazer o procedimento "na mão" (não é trivial corrigir os erros-padrão). Na especficiação abaixo, a primeira equação modela o processo de seleção. A segunda é o modelo populacional de interesse. Ambas as estimações diferem apenas na opção do "método": o de máxima verossimilhança produz erros-padrão menores (assintoticamente).
modeloHeckman_2S = selection(
  inlf ~ educ + exper + nwifeinc + age + kidslt6 + kidsge6,
  log(wage) ~ educ + exper + I(exper^2),
  data = mroz,
  method = '2step'
)
summary(modeloHeckman_2S)

modeloHeckman_ML = selection(
  inlf ~ educ + exper + nwifeinc + age + kidslt6 + kidsge6,
  log(wage) ~ educ + exper + I(exper^2),
  data = mroz,
  method = 'ml'
)
summary(modeloHeckman_ML)

# Tabela ilustrando a "igualdade" das estimativas. (Desvios decorrem, certamente, dos algoritmos de solução dos problemas que definem os respectivos métodos de estimação.)
cbind(
  '2S' = coef(modeloHeckman_2S)[-12:-14],
  'ML' = coef(modeloHeckman_ML)[-12:-13],
  'Diferença' = coef(modeloHeckman_2S)[-12:-14] - coef(modeloHeckman_ML)[-12:-13]
)
