rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)

# Abre a base de dados "crime1" do Wooldridge.
# J. Grogger (1991), “Certainty vs. Severity of Punishment,” Economic Inquiry 29, 297-309.
data(crime1)

# narr86: # times arrested, 1986
# nfarr86: # felony arrests, 1986
# nparr86: # property crme arr., 1986
# pcnv: proportion of prior convictions
# avgsen: avg sentence length, mos.
# tottime: time in prison since 18 (mos.)
# ptime86: mos. in prison during 1986
# qemp86: # quarters employed, 1986
# inc86: legal income, 1986, $100s
# durat: recent unemp duration
# black: =1 if black
# hispan: =1 if Hispanic
# born60: =1 if born in 1960
# pcnvsq: pcnv^2
# pt86sq: ptime86^2
# inc86sq: inc86^2

### SUMÁRIO ###
# 1 -- ESTIMAÇÃO
# 3 -- MEDIDA DE AJUSTE (R^2)
# 4 -- INFERÊNCIA

# estimar lpm, estimar exponencial, efeito do exp, 

# ----------------------------------------------------------------------
# POISSON -- ESTIMAÇÃO
# ----------------------------------------------------------------------
# Gráfico motivacional: é evidente a natureza discreta e limitada a poucos valores da variável de interesse.
ggplot(crime1, aes(x = narr86)) +
  geom_histogram(colour = 'black', fill = 'white', binwidth = 1) +
  geom_density(alpha = 0.1, fill= 'red') + 
  ylab('Contagem')

sort(unique(crime1$narr86))

# Modelo linear. Lembre de avaliar a significância estatística dos coeficientes considerando a heterocedasticidade dos resíduos.
modeloLinear = lm(
  narr86 ~ pcnv + avgsen + tottime + ptime86 + qemp86 + inc86 + black + hispan + born60,
  data = crime1
)

coeftest(modeloLinear, vcov. = hccm)

# Estima o modelo Poisson.
modeloPoisson = glm(
  narr86 ~ pcnv + avgsen + tottime + ptime86 + qemp86 + inc86 + black + hispan + born60,
  data = crime1,
  family = poisson
)

# Estima o modelo Quasi-Poisson (veja referência a ele no Wooldrige; o ponto é permitir alguma flexibilidade à variância da variável dependente).
modeloQuasiPoisson = glm(
  narr86 ~ pcnv + avgsen + tottime + ptime86 + qemp86 + inc86 + black + hispan + born60,
  data = crime1,
  family = quasipoisson
)

# Observe que a diferença entre os modelos é, como esperado, inexistente do lado das estimativas para os coeficientes de interesse. A distinção só aparece nos erros-padrão.
summaryPoisson = summary(modeloPoisson)
summaryQuasiPoisson = summary(modeloQuasiPoisson)

cbind(
  'Poisson' = coef(summaryPoisson)[,2],
  'Quasi-Poisson' = coef(summaryQuasiPoisson)[,2]
)

# Mais, observe que a diferença entre ambos é uma constante.
coef(summaryQuasiPoisson)[,2] / coef(summaryPoisson)[,2]

# Esta constante (sigma) pode ser estimada consistentemente "na mão" com facilidade a partir do modelo Poisson "comum" (isto é, na prática, não precisaríamos que o 'glm' fizesse o trabalho por nós).
yHat = modeloPoisson[['fitted.values']]
u2 = (crime1$narr86 - yHat)^2
sigma2 = (modeloPoisson[['df.residual']])^(-1) * sum(u2 / yHat)
sqrt(sigma2)

cbind(
  'Poisson' = sqrt(sigma2) * coef(summaryPoisson)[,2],
  'Quasi-Poisson' = coef(summaryQuasiPoisson)[,2]
)

# Lembre-se de que a interpretação dos coeficientes do modelo Poisson é simétrica ao caso de se ter um modelo log-nível (100 vezes a estimativa aproxima o efeito percentual sobre a esperança condicional da variável dependente para um incremento unitário na respectiva variável independente).
cbind(
  'Linear' = coef(modeloLinear),
  '(Quasi-)Poisson' = coef(modeloPoisson)
)

# ----------------------------------------------------------------------
# POISSON -- MEDIDA DE AJUSTE (R^2)
# ----------------------------------------------------------------------
# Aqui, o R^2 é análogo àquele do script do Tobit.
R2Poisson = cor(crime1$narr86, yHat)^2
R2Poisson

# Gráficos ilustrando o (péssimo) ajuste do modelo.
library(gridExtra)
  plot1 = ggplot(crime1, aes(x = narr86)) +
    geom_histogram(colour = 'black', fill = 'white', bins = 10) +
    ylab('Contagem')
  
  plot2 = ggplot(mapping = aes(yHat)) +
    geom_histogram(colour = 'black', fill = 'white', bins = 10) +
    xlab('Fit') +
    ylab('Contagem')
    
  grid.arrange(plot1, plot2, ncol = 2)

# ----------------------------------------------------------------------
# POISSON -- INFERÊNCIA
# ----------------------------------------------------------------------
# Na medida em que modelos Poisson são estimados por máxima verossimilhança, valem todos os comentários e procedimentos para inferência discutidos para os modelos Probit e Logit (podem usar os mesmos comandos, automáticos ou manuais).
  
# A única diferença é que vocês vão encontrar (naturalmente) diferenças se usarem os modelos Poisson ou Quasi-Poisson. O Wooldridge indica um modo de mapear a estatística LR do Poisson a do Quasi-Poisson (bastaria dividir a primeira por 'sqrt(sigma2)' calculado acima). De qualquer modo, aqui é referível usar os comandos "prontos" com os "modelos certos".
