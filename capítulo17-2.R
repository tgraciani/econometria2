rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)
library(AER) # Atenção!
library(censReg) # Atenção!

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
# 1 -- ESTIMAÇÃO E EFEITOS MARGINAIS (PEA AUTOMATIZADO)
# 2 -- EFEITOS MARGINAIS (PEA AUTOMATIZADO)
# 3 -- EFEITOS MARGINAIS (PEA MANUAL)
# 4 -- EFEITOS MARGINAIS (APE MANUAL)
# 5 -- MEDIDA DE AJUSTE (R^2)
# 6 -- INFERÊNCIA

# ----------------------------------------------------------------------
# TOBIT -- ESTIMAÇÃO
# ----------------------------------------------------------------------
# Gráfico motivacional: é evidente o "ponto de massa" na densidade da variável de horas de trabalho.
ggplot(mroz, aes(x = hours)) +
  geom_histogram(aes(y = ..density..), colour = 'black', fill = 'white', bins = 30) +
  geom_density(alpha = 0.1, fill= 'red') + 
  ylab('Densidade')

# Contagem dos casos em que 'hours' é 0 ou assume valor positivo.
mroz %>% mutate(I = if_else(hours == 0, 'Zero', 'Positivo')) %>% count(I)

# Estimação do modelo Tobit.
modeloTobit = tobit(
  hours ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6,
  data = mroz
)

# Lembre-se de que as estimativas aqui não são de interpretação imediata: elas mensuram efeitos marginais sobre a esperança condicional da variável latente (não é nosso interesse).
# OBSERVAÇÃO: Apesar de o summary retornar uma suposta estatística "t", esta é uma "z". (Calcule o p-valor na mão para veriricar.) Uma opção "correta" é o comando 'tobit' do pacote AER, mas ele não produz efeitos marginais automáticos, por isso preferi o 
summary(modeloTobit)

# Estimação de um modelo linear (apenas para comparação).
modeloLinear = lm(
  hours ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6,
  data = mroz
)

# Para analisar erros-padrão (em exercícios de inferência), é preciso levar em consideração da heterocedasticidade natural implicada pelos dados (não é possível que a distribuição condicional da variável dependente seja normal quando temos um "ponto de acumulação").
coeftest(modeloLinear, vcov. = hccm)

# Observem que o LM, como na "seção" anterior, também pode gerar previsões absurdas (no caso, negativas).
mroz$LM = modeloLinear[["fitted.values"]]

ggplot(mroz, aes(x = LM)) +
  geom_histogram(colour = 'black', fill = 'white', bins = 30) +
  xlab('Previsão / LM') + 
  ylab('Contagem')

# ----------------------------------------------------------------------
# TOBIT -- EFEITOS MARGINAIS (PEA AUTOMATIZADO)
# ----------------------------------------------------------------------
# Como a comparação das estimativas dos modelos Linear / Tobit não são diretas (são efeitos marginais sobre a esperança condicional da variável latente), é preciso fazer ajustes: no caso, "retificar" as do Tobit por um fator de correção. A alternativa mais simples é o PEA, definido como nos modelos Probit / Logit. Sua implementação é automática.
# OBSERVAÇÃO: Minha única "reclamação" com o censReg é ele reportar os testes de hipótese como "t", não "z". Mas todos os números estão corretos, tanto faz usar este pacote ou o AER.
modeloTobit_v2 = censReg(
  hours ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6,
  data = mroz
)
summary(modeloTobit_v2)

margEff(modeloTobit_v2)

# Tabelas comparativas. Estas poderiam ser "melhor automatizadas" trocando o "-8" por "-length(coef(modeloTobit))".
cbind('Modelo Linear' = coef(modeloLinear), 'Tobit' = coef(modeloTobit_v2)[-8])

cbind('Modelo Linear' = coef(modeloLinear)[-1], 'Tobit' = margEff(modeloTobit_v2))

# ----------------------------------------------------------------------
# TOBIT -- EFEITOS MARGINAIS (PEA MANUAL)
# ----------------------------------------------------------------------
# Salva informações úteis: os coeficientes do Tobit em um vetor (mas excluindo o sigma) e o sigma, mas corridigo (o comando de estimação retorna o logaritmo de sigma, então é preciso reverter a operação).
B = coef(modeloTobit)
sigma = modeloTobit[['scale']]

# Cria o indivíduo médio baseado nas observações da amostra.
iMédio = mroz %>% summarise_at(
  vars(nwifeinc, educ, exper, age, kidslt6, kidsge6), 
  ~ mean(.x, na.rm = TRUE)
)
iMédio

# Calcula o fator de correção do PEA.
fatorPEA = as.numeric(
  pnorm(
    (B[1] + sum(B[2:7] * iMédio)) / sigma
  )
)
fatorPEA

# Calcula os efeitos marginais e checagem que batem com o comando automatizado. (Vejam que dá uma diferença na sexta casa decimal da última variável; deve ser de uma diferença minúscula entre as estimativas dos dois comandos, que usam métodos distintos para resolver o proxima de maximizar a log-verossimilhança.)
efeitosPEA = B * fatorPEA

cbind('Tobit Automático' = margEff(modeloTobit_v2), 'Tobit Manual' = efeitosPEA[-1])

# ----------------------------------------------------------------------
# TOBIT -- CÁLCULO MANUAL DO APE
# ----------------------------------------------------------------------
# Adiciona alguns cálculos essenciais à base de dados: (i) o "XB" de cada observação e (ii) o valor da função de distribuição da normal avaliada nesse "XB".
# OBSERVAÇÃO: Veja que o comando 'tobit' já calcula o 'XB' (abra o objeto resultado do modelo e veja os 'linear.predictors'). Preferi fazer na mão caso alguém decida ficar com o 'censReg', que não faz essa conta sozinho.
mroz = mroz %>%
  mutate(
    XB = B[1] + (B[2] * nwifeinc) + (B[3] * educ) + (B[4] * exper) + (B[5] * age) + (B[6] * kidslt6) + (B[7] * kidsge6),
    
    fatorIndividual = pnorm(XB / sigma)
  )

# Para obter o fator do APE, basta tomar a média de (ii) calculado acima.
fatorAPE = mean(mroz$fatorIndividual)
fatorAPE

# Calcula os efeitos marginais e faz a comparação com resultados anteriores.
efeitosAPE = B * fatorAPE

cbind('Modelo Linear' = coef(modeloLinear)[-1],
  'Tobit (PEA)' = efeitosPEA[-1],
  'Tobit (APE)' = efeitosAPE[-1])

# ----------------------------------------------------------------------
# TOBIT -- MEDIDA DE AJUSTE (R^2)
# ----------------------------------------------------------------------
# Para calcular a medida de R^2 sugerida pelo Wooldridge (correlação entre a variável dependente e o fit do modelo), é preciso calcular a esperança condicional da variável de interesse (que é função não linear dos regressores e dos coeficientes do modelo).
mroz = mroz %>% mutate(
  yHat = (pnorm(XB / sigma) * XB)  + (sigma * (dnorm(XB / sigma)))
)

R2Tobit = cor(mroz$hours, mroz$yHat)^2
R2Tobit

# SComparação com o modelo linear
R2Linear = cor(mroz$hours, modeloLinear[["fitted.values"]])^2
R2Linear

# Só para visualizar as previsões do modelo e ilustrar que, diferentemente do LM, essas são não negativas.
ggplot(mroz, aes(x = yHat)) +
  geom_histogram(colour = 'black', fill = 'white', bins = 30) +
  xlab('Previsão / Tobit') + 
  ylab('Contagem')

# ----------------------------------------------------------------------
# TOBIT -- INFERÊNCIA
# ----------------------------------------------------------------------
# Na medida em que modelos Tobit são estimados por máxima verossimilhança, valem todos os comentários e procedimentos para inferência discutidos para os modelos Probit e Logit (podem usar os mesmos comandos, automáticos ou manuais).
