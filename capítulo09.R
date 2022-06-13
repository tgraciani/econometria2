rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)

### SUMÁRIO ###
# 1 -- TESTE RESET PARA ESPECIFICAÇÃO FUNCIONAL
# 2 -- TESTES PARA ALTERNATIVAS NÃO ANINHADAS
# 3 -- SIMULAÇÃO: ERRO DE MEDIDA NA VARIÁVEL DEPENDENTE
# 4 -- SIMULAÇÃO: ERRO DE MEDIDA NA VARIÁVEL INDEPENDENTE

# Consulte <https://cran.r-project.org/web/packages/wooldridge/wooldridge.pdf> para obter um dicionário das bases de dados do Wooldridge.

# ----------------------------------------------------------------------
# TESTE RESET PARA ESPECIFICAÇÃO FUNCIONAL
# ----------------------------------------------------------------------
data(hprice1)

# price: house price, $1000s
# lotsize: size of lot in square feet
# sqrft: size of house in square feet
# bdrms: number of bdrms

# Implementação manual do teste RESET. Basta (i) adicionar à regressão um polinômio do "fit" do modelo populacional e fazer um teste de exclusão conjunta.
mod1 = lm(
  price ~ lotsize + sqrft + bdrms,
  data = hprice1
)

hprice1$yHat1 = fitted(mod1)

modReset = lm(
  price ~ lotsize + sqrft + bdrms + I(yHat1^2) + I(yHat1^3),
  data = hprice1
)

linearHypothesis(modReset, c('I(yHat1^2)', 'I(yHat1^3)'))

# Observe a correta distribuição da estatística de teste (é necessário considerar a "adição de regressores").
pf(4.6682, 2, 88-4-2, lower.tail = FALSE)

# Implementação automática do teste RESET, do pacote 'lmtest'.
resettest(mod1, power = 2:3)

# ----------------------------------------------------------------------
# TESTES PARA ALTERNATIVAS NÃO-ANINHADAS
# ----------------------------------------------------------------------
# Teste para alternativas não aninhadas (de Mizon-Richard). É necessário apenas estimar um "modelo conjunto" e fazer um teste de exclusão (conjunta). Observe os dois "lados" possíveis do teste.
mod2 = lm(
  price ~ log(lotsize) + log(sqrft) + bdrms,
  data = hprice1
)

summary(mod1)
summary(mod2)

modMR = lm(
  price ~ lotsize + sqrft + bdrms + log(lotsize) + log(sqrft),
  data = hprice1
)
summary(modMR)

linearHypothesis(modMR, c('log(lotsize)', 'log(sqrft)'))
linearHypothesis(modMR, c('lotsize', 'sqrft'))

# Condução automatizada do teste, habilitada pelo pacote 'lmtest'.
encomptest(mod1, mod2, data = hprice1)

# Cálculo te outro teste, agora com a abordagem de Davidson-MacKinnon. Similarmente ao teste RESET, basta adicionar o "fit" do modelo alternativo e conduzir um teste de exclusão. 
hprice1$yHat2 = fitted(mod2)

modDM = lm(
  price ~ lotsize + sqrft + bdrms + yHat2,
  data = hprice1
)
coeftest(modDM)

modDM = lm(
  price ~ log(lotsize) + log(sqrft) + bdrms + yHat1,
  data = hprice1
)
coeftest(modDM)

# ----------------------------------------------------------------------
# EXEMPLO DE "FOR LOOP"
# ----------------------------------------------------------------------
# Criação básica de sequências.
1:5
seq(from = 0, to = 100, by = 7)

# Exemplo de loop.
for (i in 1:5) {
  print(paste('Neste passo da contagem, o valor de i é ', i, '.', sep = ''), quote = FALSE)
}

# ----------------------------------------------------------------------
# SIMULAÇÃO: ERRO DE MEDIDA NA VARIÁVEL DEPENDENTE
# ----------------------------------------------------------------------
rm(list = ls())

# Escolha dos valores verdadeiros dos parâmetros do modelo populacional de interesse.
B0 = 5
B1 = -2.5

# Criação de quatro vetores "vazios". Cada par carregará, respectivamente, estimativas de OLS sem e com erro de medida.
B0Hat = numeric(1000)
B1Hat = numeric(1000)
B0HatME = numeric(1000)
B1HatME = numeric(1000)

# Criação da nossa variável independente (extrações de uma distribuição normal de média e erro padrão 3).
x = rnorm(1000, mean = 3, sd = 3)
mean(x)

# Primeiro, no loop se criam erros e variável dependente para, na sequência, estimar o modelo populacional por OLS, guardando a estimativa de interesse. Segundo, refaz-se o exercício, mas adicionando um erro de medida 'e0' sobre a variável dependente.
for (i in 1:1000) {
  u = rnorm(1000)
  yStar = B0 + B1*x + u
  BHat = coef(lm(yStar ~ x))
  B1Hat[i] = BHat['x']
  
  e0 = rnorm(1000)
  y = yStar + e0
  BHatME = coef(lm(y ~ x))
  B1HatME[i] = BHatME['x']
}

# Comparação das estimativas (do coeficiente de 'x').
cbind(
  'Com Erro {A}' = mean(B1HatME),
  'Sem Erro {B}' = mean(B1Hat),
  '{A}/{B}' = mean(B1HatME)/mean(B1Hat)
)

# Comparação das variâncias das estimativas.
cbind(
  'Com Erro {A}' = var(B1HatME),
  'Sem Erro {B}' = var(B1Hat),
  '{A}/{B}' = var(B1HatME)/var(B1Hat)
)

# Refaz o exercício de simulação, mas agora considerando um erro de medida cuja média não é zero.
for (i in 1:1000) {
  u = rnorm(1000)
  yStar = B0 + B1*x + u
  BHat = coef(lm(yStar ~ x))
  B0Hat[i] = BHat['(Intercept)']
  B1Hat[i] = BHat['x']
  
  e0 = rnorm(1000, mean = 7)
  y = yStar + e0
  BHatME = coef(lm(y ~ x))
  B0HatME[i] = BHatME['(Intercept)']
  B1HatME[i] = BHatME['x']
}

cbind(
  'Com Erro {A}' = mean(B1HatME),
  'Sem Erro {B}' = mean(B1Hat),
  '{A}/{B}' = mean(B1HatME)/mean(B1Hat)
)

cbind(
  'Com Erro {A}' = var(B1HatME),
  'Sem Erro {B}' = var(B1Hat),
  '{A}/{B}' = var(B1HatME)/var(B1Hat)
)

cbind(
  'Com Erro {A}' = mean(B0HatME),
  'Sem Erro {B}' = mean(B0Hat),
  '{A}-{B}' = mean(B0HatME) - mean(B0Hat)
)

# Refaz o primeiro exercício de simulação, mas considerando um erro de medida minimamente correlacionado com o regressor. (Experimente variar essa correlação)
for (i in 1:1000) {
  u = rnorm(1000)
  yStar = B0 + B1*x + u
  BHat = coef(lm(yStar ~ x))
  B0Hat[i] = BHat['(Intercept)']
  B1Hat[i] = BHat['x']
  
  e0 = rnorm(1000) + (0.1 * x)
  y = yStar + e0
  BHatME = coef(lm(y ~ x))
  B0HatME[i] = BHatME['(Intercept)']
  B1HatME[i] = BHatME['x']
}

cbind(
  'Com Erro {A}' = mean(B1HatME),
  'Sem Erro {B}' = mean(B1Hat),
  '{A}/{B}' = mean(B1HatME)/mean(B1Hat)
)

cbind(
  'Com Erro {A}' = var(B1HatME),
  'Sem Erro {B}' = var(B1Hat),
  '{A}/{B}' = var(B1HatME)/var(B1Hat)
)

cbind(
  'Com Erro {A}' = mean(B0HatME),
  'Sem Erro {B}' = mean(B0Hat),
  '{A}/{B}' = mean(B0HatME)/mean(B0Hat)
)

# ----------------------------------------------------------------------
# SIMULAÇÃO: ERRO DE MEDIDA NA VARIÁVEL INDEPENDENTE
# ----------------------------------------------------------------------
rm(list = ls())

# Define parâmetros e prepara vetores para o exercício de simulação.
B0 = 5
B1 = -2.5

B1Hat = numeric(1000)
B1HatME = numeric(1000)

# Aqui, o exercício segue a mesma lógica do anterior. A diferença relevante é a observação do regressor, 'x'. No caso, a simulação cria 'xStar' e 'e0' independentemente, de modo que, por hipótese, 'e0' é correlacionado com 'x', nosso regressor com erro de medida (é a hipótese de "Classical Error-in-Variables", CEV).
# Refaça o exercício aumentando a variância de 'xStar' para chegar à conclusão de que, assintoticamente, o viés de atenuação é tanto menor quanto maior for a variância de 'xStar' relativamente àquela do erro.
#xStar = rnorm(1000, mean = 3, sd = 3)
#xStar = rnorm(1000, mean = 3, sd = 3*2)
xStar = rnorm(1000, mean = 3, sd = 3/2)

for (i in 1:1000) {
  u = rnorm(1000)
  y = B0 + B1*xStar + u
  BHat = coef(lm(y ~ xStar))
  B1Hat[i] = BHat['xStar']
  
  e0 = rnorm(1000)
  x = xStar + e0
  BHatME = coef(lm(y ~ x))
  B1HatME[i] = BHatME['x']
}

cbind(
  'Com Erro {A}' = mean(B1HatME),
  'Sem Erro {B}' = mean(B1Hat),
  '{A}/{B}' = mean(B1HatME)/mean(B1Hat)
)

cbind(
  'Com Erro {A}' = var(B1HatME),
  'Sem Erro {B}' = var(B1Hat),
  '{A}/{B}' = var(B1HatME)/var(B1Hat)
)

# Refaz o exercício, mas agora considerando o caso em que 'x' (com erro) é não correlacionado com 'e0', de modo que o valor sem erro da variável é correlacionado com o problema de medida.
x = rnorm(1000, mean = 3, sd = 3)

for (i in 1:1000) {
  u = rnorm(1000)
  e0 = rnorm(1000)
  xStar = x - e0
  y = B0 + B1*xStar + u
  BHat = coef(lm(y ~ xStar))
  B1Hat[i] = BHat['xStar']
  
  BHatME = coef(lm(y ~ x))
  B1HatME[i] = BHatME['x']
}

cbind(
  'Com Erro {A}' = mean(B1HatME),
  'Sem Erro {B}' = mean(B1Hat),
  '{A}/{B}' = mean(B1HatME)/mean(B1Hat)
)

cbind(
  'Com Erro {A}' = var(B1HatME),
  'Sem Erro {B}' = var(B1Hat),
  '{A}/{B}' = var(B1HatME)/var(B1Hat)
)
