rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)

# Abre a base de dados "mroz" do Wooldridge.
# T. A. Mroz (1987), “The Sensitivity of an Empirical Model of Married Women’s Hours of Work to Economic and Statistical Assumptions,” Econometrica 55, 765-799.
data(mroz)

# [...] husband’s earnings (nwifeinc, measured in thousands of dollars), years of education (educ), past years of labor market experience (exper), age, number of children less than six years old (kidslt6), and number of kids between 6 and 18 years of age (kidsge6).



# ----------------------------------------------------------------------
# LPM -- ESTIMAÇÃO E INFERÊNCIA
# ----------------------------------------------------------------------
# Cria uma variável igual ao quadrado dos anos de experiência. (Só para ilustrar novamente a criação de variáveis, não vamos usar nos modelos dos exemplos.)
mroz = mroz %>% mutate(exper2 = exper ^ 2)
mroz = mutate(mroz, exper2 = exper ^2)
mroz$exper2 = mroz$exper ^ 2

# Estima um LPM usando (i) uma fórmula e (ii) especificando uma base de dados.
modeloLPM = lm(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz)

# Exibe os resultados do modelo LPM.
summary(modeloLPM)

# Testes de hipótese individuais com estimação robusta dos erros padrão.
coeftest(modeloLPM, vcov. = hccm)

# Cria duas "observações hipotéticas" para a qual se pode calcular as probabilidades previstas pelo LPM. O objetivo é mostrar previsões absurdas do modelo.
mulheresRadicais = list(
  nwifeinc = c(100, 0), #Salário do marido de uma é 100 e de outra, 0.
  educ = c(5, 17), # Anos de educação de uma é 5 e de outra, 17.
  exper = c(0, 30), # Assim sucessivamente...
  age = c(20, 52),
  kidslt6 = c(2, 0),
  kidsge6 = c(0, 0)
)

# Calcula a probabilidade prevista usando as observações hipotéticas com os resultados do modelo LPM estimado.
predict(modeloLPM, mulheresRadicais)

# Para só uma observação, pode ser mais simples usar o formato abaixo. Na primeira linha "salvamos" os coeficientes. Na segunda, dizemos quais são os valores das variáveis dependentes para nossa observação fictícia. Na terceira, multiplicamos as variáveis pelos respectivos coeficientes, somando os resultados para chegar na probabilidade estimada.
coeficientesLPM = modeloLPM[["coefficients"]]
mulherRadical1 = c(1, 100, 5, 0, 20, 2 , 0)
sum(mulherRadical1 * coeficientesLPM)

# Cria duas colunas no banco de dados. Uma é apenas uma etiqueta das observações (para fazer gráficos abaixo). A outra salva os valores ajustados pelo modelo na própria base de dados (como se fosse uma variável).
mroz$ID = as.numeric(row.names(mroz))
mroz$fit = modeloLPM[["fitted.values"]]

# ----------------------------------------------------------------------
# LPM -- MEDIDAS DE AJUSTE
# ----------------------------------------------------------------------
# Cria três variáveis. Primeiro, uma dummy indicando se o LPM prevê que a o valor verdadeiro da variável dependente seja 1. Segundo, só comparamos a primeira variável com a dependente observada, chechando se o modelo fez uma previsão correta ou não. Terceiro, indicamos quais observações produziram probabilidades absurdas.
mroz = mroz %>% mutate(
  resultadoPrevisto = if_else(fit >= 0.5, 1, 0),
  indicadorAcerto = if_else(resultadoPrevisto == inlf, 'Correto', 'Errado'),
  indicadorAbsurdo = if_else(fit < 0 | fit > 1, 'Problemático', 'OK')
)

# Proporção de erros e acertos (comparando previsões e dados observados).
mroz %>% count(indicadorAcerto) %>% mutate(freq = n / sum(n))

# Proporção probabilidades absurdas.
mroz %>% count(indicadorAbsurdo) %>% mutate(freq = n / sum(n))

# ----------------------------------------------------------------------
# LPM -- GRÁFICOS (OPCIONAL)
# ----------------------------------------------------------------------
# Gráfico ilustrando acertos e erros.
ggplot(mroz, aes(x = age, y = fit, color = indicadorAcerto)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.5) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black', size = 0.5) +
  ylab('Probabilidade Ajustada')

# Gráfico ilustrando a "distribuição" que o modelo gera.
ggplot(mroz, aes(x = fit)) +
  geom_histogram(color = 'black', fill = 'white', binwidth = 0.05) +
  xlab('Probabilidade Ajustada') +
  ylab('Contagem')

# Versões mais simples/diretas (mas menos customizáveis) dos gráficos acima.
plot(mroz$fit)
hist(mroz$fit)

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- DENSIDADES (OPCIONAL)
# ----------------------------------------------------------------------
# Calcula o valor da densidade para argumentos quaisquer.
dnorm(2)
dlogis(1.23)

# Cria um banco de dados vazio com 10000 linhas na primeira linha. Na segunda, subdivide o intervalo [-5, 5] em 10000 partes equidistantes. Baseados nestes valores, as últimas linhas calculam o valor das densidades das distribuições normal padrão e logística.
dd = data.frame(matrix(nrow = 10000, ncol = 0))
dd$sequência = seq(-5, 5, length = 10000)
dd$distribuiçãoNormal = dnorm(dd$sequência)
dd$distribuiçãoLogística = dlogis(dd$sequência)

# Faz o gráfico com as 10000 "observações" das distribuições.
ggplot(data = dd, aes(x = sequência)) +
  geom_line(aes(y = distribuiçãoNormal), color = 'red') +
  geom_line(aes(y = distribuiçãoLogística), color = 'blue') +
  ylab('Densidade') +
  xlab('')

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- ESTIMAÇÃO
# ----------------------------------------------------------------------
# Estima modelos probit e logit. Observação importante: aqui, erros-padrão e testes de hipótese individuais já são assintóticos.
modeloProbit = glm(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, family = binomial(link=probit))
summary(modeloProbit)

modeloLogit = glm(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, family = binomial(link=logit))
summary(modeloLogit)

# Monta uma tabelinha comparando as estimativas dos dois modelos.
cbind(LPM = modeloLPM[["coefficients"]], Probit = modeloProbit[["coefficients"]], Logit = modeloLogit[["coefficients"]])

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- INFERÊNCIA ASSINTÓTICA: TESTE LR
# ----------------------------------------------------------------------
# Aqui, vou discutir as coisas para o modelo probit, mas tudo poderia ser feito para o logit só mudando os argumentos das funções.
# Para testes de exclusão múltipla, uma alternativa é o teste LR. Nas duas linhas abaixo, se testa o modelo contra a exclusão de todas as covariadas.
lrtest(modeloProbit)

# Estima um novo modelo probit, mas sem as variáveis de filhos. É o nosso "modelo restrito". 
modeloProbit = glm(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, family = binomial(link=probit))
modeloProbitRestrito = glm(inlf ~ nwifeinc + educ + exper + age, data = mroz, family = binomial(link=probit))

# Faz o teste de exclusão conjunta das variáveis de filhos comparando os dois modelos, irrestrito e restrito. (A hipótese nula é de que os modelos são equivalentes.)
lrtest(modeloProbit, modeloProbitRestrito)

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- INFERÊNCIA ASSINTÓTICA: TESTE WALD
# ----------------------------------------------------------------------
# A implementação é direta com a biblioteca "lmtest":
waldtest(modeloProbit, c('kidslt6', 'kidsge6'), test = c("Chisq"))

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- INFERÊNCIA ASSINTÓTICA: TESTE LM
# ----------------------------------------------------------------------
# Vamos implementar uma "versão esperta" do teste LM (não achei pacotes que façam isso sozinhos; veja a p. 570 do Wooldridge de pós-graduação, se tiver curiosidade). Para tanto, é necessário seguir os seguintes passos:
# (1) Estime o modelo restrito (já feito acima)
summary(modeloProbitRestrito)

# (2) Salve os resíduos da estimação. Também salve os valores de "G" e "g" (probabilidades previstas e "preditores lineares") para cada observação. Vou fazer isso num novo dataframe para evitar confusões.
dd = mroz
dd = dd %>% mutate(
  u = modeloProbitRestrito[["residuals"]],
  G = modeloProbitRestrito[["fitted.values"]],
  g = modeloProbitRestrito[["linear.predictors"]],
  correção_u = 1 / sqrt(G * (1-G)),
  correção = g / sqrt(G * (1-G))
)

# (3) Transforme as variáveis de acordo com o descrito no Wooldridge de pós-graduação.
dd = dd %>% transmute(
  u = u * correção_u,
  nwifeinc = nwifeinc * correção,
  educ = educ * correção,
  exper = exper * correção,
  age = age * correção,
  kidslt6 = kidslt6 * correção,
  kidsge6 = kidsge6 * correção
)

# (4) Regrida os resíduos no vetor completo de regressores, com todas as variáveis (resíduos inclusive) ponderadas por OLS.
OLSauxiliar = lm(u ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = dd)

# (5) Calcule o R2 (descentralizado) do modelo. Multiplique pelo tamanho da amostra. Calcule o p-valor associado.
library(rsq)
LM = nobs(OLSauxiliar) * rsq(OLSauxiliar, adj = FALSE)
LM
pchisq(LM, 2, lower.tail = FALSE)

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- MEDIDAS DE AJUSTE
# ----------------------------------------------------------------------
# Por padrão, a log-verossimilhança da estimação não é reportada, mas poder ser recuperada com a função "logLik".
logLik(modeloProbit)
logLik(modeloProbitRestrito)
logLik(modeloLogit)

# É útil obter o valor dessa log-verossimilhança para calcular o pseudo R2 de McFaden.
1 - (modeloProbit$deviance / modeloProbit$null.deviance)
1 - (modeloLogit$deviance / modeloLogit$null.deviance)

# Salva na base de dados o valor ajustado pelo modelo, a previsão feita para a variável dependente e um indicador para o caso de a previsão estar correta. Essas novas colunas são criadas tanto para o modelo probit como o logit.
mroz = mroz %>% mutate(
  fitProbit = modeloProbit[["fitted.values"]],
  previsãoProbit = if_else(fitProbit >= 0.5, 1, 0),
  acertosProbit = if_else(previsãoProbit == inlf, 1, 0),
  
  fitLogit = modeloLogit[["fitted.values"]],
  previsãoLogit = if_else(fitLogit >= 0.5, 1, 0),
  acertosLogit = if_else(previsãoLogit == inlf, 1, 0)
)

# Comparação do desempenho dos dois modelos.
mroz %>% count(previsãoProbit, previsãoLogit) %>% mutate(freq = n / sum(n))

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- EFEITOS MARGINAIS
# ----------------------------------------------------------------------
# "Salva" os coeficientes do modelo probit em um novo objeto
coeficientesProbit = coef(modeloProbit)

# Define um indivíduo alternativo (indicando os valores das variáveis dependentes).
indivíduoAntes = c(1, 20, 6, 10, 42, 0, 0)
indivíduoDepois = c(1, 20, 6, 10, 42, 1, 0)

# Calcula o efeito marginal desse indivíduo (do "modo aproximado").
dnorm(sum(indivíduoAntes * coeficientesProbit)) * coeficientesProbit

# Calcula o efeito marginal desse indivíduo (do "modo exato").
pnorm(sum(indivíduoDepois * coeficientesProbit)) - pnorm(sum(indivíduoAntes * coeficientesProbit))

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- PEA e APE
# ----------------------------------------------------------------------
# Determina o "indivíduo médio" da amostra.
indivíduoMédio = mroz %>% summarise_at(
  vars(nwifeinc, educ, exper, age, kidslt6, kidsge6), 
  ~ mean(.x, na.rm = TRUE)
)
indivíduoMédio

# Calcula o PEA.
peaProbit = coef(modeloProbit) * dnorm(predict(modeloProbit, indivíduoMédio))
peaLogit = coef(modeloLogit) * dlogis(predict(modeloLogit, indivíduoMédio))

# Exibe tabela comparativa.
cbind('PEA Probit' = peaProbit, 'PEA Logit' = peaLogit)

# Salva os valores de Xβ do modelo em objetos à parte.
xbProbit = predict(modeloProbit)
xbLogit = predict(modeloLogit)

# Calcula o APE.
apeProbit = coef(modeloProbit) * mean(dnorm(xbProbit))
apeProbit = mean(coef(modeloProbit) * dnorm(xbProbit))
apeLogit = coef(modeloLogit) * mean(dlogis(xbLogit))

# Exibe tabela comparativa.
cbind('APE Probit' = apeProbit, 'APE Logit' = apeLogit)

# ----------------------------------------------------------------------
# PROBIT E LOGIT -- PEA/APE NO PACOTE 'MFX' (OPCIONAL)
# ----------------------------------------------------------------------
# Se não tiver instalado o pacote, execute a linha abaixo.
install.packages('mfx')

# Carrega a biblioteca.
library(mfx)

# Cálculo automático do PEA.
probitmfx(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, atmean = TRUE)
logitmfx(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, atmean = TRUE)

# Cálculo automático do APE.
probitmfx(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, atmean = FALSE)
logitmfx(inlf ~ nwifeinc + educ + exper + age + kidslt6 + kidsge6, data = mroz, atmean = FALSE)