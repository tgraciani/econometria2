rm(list = ls())
#options(scipen = 100, digits = 7)
source('~/Dropbox/doutorado/monitoria/monitorias/display.R')

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)
library(AER) # Atenção!

### SUMÁRIO ###
# 1 -- MÉTODO DE INSTRUMENTAL VARIABLES (IV)
# 3 -- MÉTODO DE TWO-SATGE LEAST SQUARES (2SLS)
# 4 -- 2SLS: ESTIMAÇÃO MANUAL I
# 5 -- 2SLS: ESTIMAÇÃO MANUAL II
# 6 -- TESTE DE ENDOGENEIDADE (DE WU-HAUSMAN)
# 7 -- TESTE DE RESTRIÇÕES SOBREIDENTIFICADAS

# Consulte <https://cran.r-project.org/web/packages/wooldridge/wooldridge.pdf> para obter um dicionário das bases de dados do Wooldridge.

# ----------------------------------------------------------------------
# MÉTODO DE INSTRUMENTAL VARIABLES (IV)
# ----------------------------------------------------------------------
data(mroz)

# educ: years of schooling
# wage: est. wage from earn, hrs
# fatheduc: father’s years of schooling
# exper: actual labor mkt exper
# lwage: log(wage)
# expersq: exper^2

# Seleciona a subamostra cuja variável 'wage' pode ser transformada pela função 'log'. 
mroz = mroz %>% filter(!is.na(lwage))

# Para introduzir o método de variáveis instrumentais, considere um modelo linear simples que explica 'lwage' por um intercepto e pela variável 'educ'. Como vocês viram, neste caso particular (de um regressor endógeno e um instrumento), as fórmulas dos estimadores de OLS e IV podem ser expressas como um quociente de covariâncias.
b_OLS = with(
  mroz,
  cov(lwage, educ) / var(educ)
)
b_OLS

b_IV = with(
  mroz,
  cov(lwage, fatheduc) / cov(educ, fatheduc)
)
b_IV

# A estimação automática do modelo por OLS já é conhecida por vocês. Para estimar via IV, existe uma implementação pronta no pacote 'AER'. Note que se indica o instrumento separando-o por uma '|' do modelo de interesse. A sintaxe do comando exige que todas as variáveis exógenas estejam presentes do lado direito de '|'.
mod_OLS = lm(
  lwage ~ educ,
  data = mroz
)
summary(mod_OLS)

mod_IV = ivreg(
  lwage ~ educ | fatheduc,
  data = mroz
)
summary(mod_IV)

# Note que, estimando o modelo usando o próprio regressor como seu instrumento, o 'ivreg' retorna o estimador de OLS.
mod_IV = ivreg(
  log(wage) ~ educ | educ,
  data = mroz
)
summary(mod_IV)

# Sob homocedasticidade, é simples calcular a variância assintótica dos estimadores de IV. Para computá-la, é necessário obter (i) um estimador para a variância condicional do erro populacional ('sigma2'), (ii) o SST de 'educ' e (iii) e o R-squared da regressão, estimada via OLS, de 'educ' contra seu instrumento, 'fatheduc'.
show_math(
  '$$\\text{Avar} \\hat{\\beta}_1^{\\text{IV}} = \\frac{\\sigma^2_{\\varepsilon}}{\\text{SST}_x \\cdot \\text{R}^2_{x,z}}$$',
  '$$\\text{Avar} \\hat{\\beta}_1^{\\text{OLS}} = \\frac{\\sigma^2_{\\varepsilon}}{\\text{SST}_x}$$'
)

sigma2 = (1 / mod_IV[['df.residual']]) * sum(mod_IV[['residuals']]^2)

SST_x = with(
  mroz,
  sum((educ - mean(educ))^2)
)

R2_xz = summary(
  lm(
    educ ~ fatheduc,
    data = mroz
  )
)[['r.squared']]

AvarIV = sigma2 / (SST_x * R2_xz)

# Que pode ser diretamente comparada com AvarOLS:
AvarOLS = (1 / mod_OLS[['df.residual']]) * sum(mod_OLS[['residuals']]^2) / SST_x

# ----------------------------------------------------------------------
# MÉTODO DE TWO-SATGE LEAST SQUARES (2SLS)
# ----------------------------------------------------------------------
rm(list = setdiff(ls(), c('show_math')))
data(mroz)

# A função 'ivreg' do 'AER' também estima modelos "mais complicados". A sintaxe é a mesma: especifique seu modelo de interesse como se fazia com o comando 'lm' e, separando por uma '|', indique a lista de variáveis exógenas (tanto regressores exógenos do modelo como as variáveis excluídas). Todas as variáveis ausentes do lado direito de '|' são interpretadas como endógenas e são, portanto, instrumentalizadas.
mod_2SLS = ivreg(
  lwage ~ educ + exper + expersq | motheduc + fatheduc + exper + expersq,
  data = mroz
)
summary(mod_2SLS)

# Com uma opção adicional, o 'summary()' já executa um teste F sobre o primeiro estágio, o teste de Wu-Hausman de endogeneidade e o de Sargan de overidentifying restrictions.
summary(mod_2SLS, diagnostics = TRUE)

# Observe que o output dos testes não muda. A correção de heterocedasticidade arbitrária é aplicada somente sobre os erros-padrão do modelo de interesse. Observe a documentação de 'vcovHC' para uma lista de opções.
summary(mod_2SLS, vcov. = vcovHC(mod_2SLS, type = 'HC3'), diagnostics = TRUE)

# ----------------------------------------------------------------------
# 2SLS: ESTIMAÇÃO MANUAL I
# ----------------------------------------------------------------------
rm(list = setdiff(ls(), c('mod_2SLS', 'show_math')))
data(mroz)

# Para obter os estimadores de 2SLS manualmente, é preciso estimar dois modelos por OLS, do primeiro e do segundo estágio. Para facilitar, primeiro restrinja a amostra àquela que, de fato, poderá ser usada nas estimações.
mroz = mroz %>% filter(!is.na(lwage))

# Então, estime o primeiro estágio por OLS.
mod_1S = lm(
  educ ~ motheduc + fatheduc + exper + expersq,
  data = mroz
)

# Na sequência, são salvos os fitted values do modelo.
mroz$educHat = mod_1S[['fitted.values']]

# Por fim, estima-se o segundo estágio substituindo a variável endógena pelos fitted values na equação de interesse.
mod_2S = lm(
  lwage ~ educHat + exper + expersq,
  data = mroz
)

# Tabela comparativa.
# OBSERVAÇÃO. Erros-padrão do segundo estágio estimado manualmente estão incorretos (é desconsiderado que, na verdade, o erro do modelo é composto).
cbind(
  'AER::ivreg' = coef(mod_2SLS),
  'Manual' = coef(mod_2S)
)

# Como estimamos o primeiro estágio manualmente, é possível reproduzir a linha "Weak instruments" do 'summary' aplicado sobre o modelo estimado com o 'ivreg'.
linearHypothesis(mod_1S, c('motheduc', 'fatheduc'))

# Talvez mais interessante ainda seja a possibilidade de conduzir o respectivo teste robusto. Observe que o Wooldridge cita Olea and Pflueger (2013) para sugerir uma regra de bolso no caso de um regressor endógeno: haveria evidência de que os instrumentos usados são "fortes" se a estatística F do teste abaixo (robusto) for igual ou maior que 20.
linearHypothesis(mod_1S, c('motheduc', 'fatheduc'), vcov. = vcovHC)

# ----------------------------------------------------------------------
# 2SLS: ESTIMAÇÃO MANUAL II
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('mod_2SLS', 'mod_1S', 'mod_2S', 'show_math')])
data(mroz)
mroz = mroz %>% filter(!is.na(lwage))

# Uma estimação manual completamente correta (incluindo erros-padrão) pode ser facilmente obtida usando as fórmulas matriciais dos estimadores de 2SLS.
show_math(
  '\\begin{bmatrix} \\beta_0 \\\\ \\beta_1 \\\\ \\vdots \\\\ \\beta_k \\end{bmatrix} = \\left( \\mathbf{X}^T \\mathbf{Z} (\\mathbf{Z}^T \\mathbf{Z})^{-1} \\mathbf{Z}^T \\mathbf{X}\\right)^{-1} \\mathbf{X}^T \\mathbf{Z} (\\mathbf{Z}^T \\mathbf{Z})^{-1} \\mathbf{Z}^T \\mathbf{y}',
  css = 'color: black; font-size:17px;'
)

# Para fazer a estimação toda "na mão", é preciso criar a variável referente ao intercepto.
mroz$Intercept = 1

# Na sequência, é preciso criar as matrizes de variáveis dependentes, variáveis do modelo e variáveis exógenas (inclusas e não inclusas no modelo).
y = data.matrix(
  mroz %>% select(lwage)
)

X = data.matrix(
  mroz %>% select(Intercept, educ, exper, expersq)
)

Z = data.matrix(
  mroz %>% select(Intercept, exper, expersq, motheduc, fatheduc)
)

# A estimação de dois estágios pode sair numa tacada só:
b_2SLS = as.vector(
  solve(t(X)%*%Z %*% solve(t(Z)%*%Z) %*% t(Z)%*%X) %*% t(X)%*%Z %*% solve(t(Z)%*%Z) %*% t(Z)%*%y
)

# Tabela comparativa.
cbind(
  'AER::ivreg' = coef(mod_2SLS),
  'Matricial' = b_2SLS
)

# Matricialmente também é fácil obter a variância assintótica dos estimadores (supondo homocedasticidade, como no output pronto do 'AER').
show_math(
  '\\text{Avar} \\hat{\\beta} = (\\mathbf{X}^T \\mathbf{Z} (\\mathbf{Z}^T \\mathbf{Z})^{-1} \\mathbf{Z}^T \\mathbf{X})^{-1} \\frac{1}{n-k} \\mathbf{e}^T \\mathbf{e}',
  css = 'color: black; font-size:20px;'
)

e = y - X %*% b_2SLS

sigma2 = (1/(dim(X)[1] - dim(X)[2])) * as.vector((t(e) %*% e))

Avar = solve(t(X)%*%Z %*% solve(t(Z)%*%Z) %*% t(Z)%*%X) * sigma2

SE_2SLS = sqrt(diag(Avar))

cbind(
  'AER::ivreg' = summary(mod_2SLS)$coefficients[,2],
  'Matricial' = SE_2SLS
)

# ----------------------------------------------------------------------
# TESTE DE ENDOGENEIDADE
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('mod_2SLS', 'mod_1S', 'mod_2S')])
data(mroz)
mroz = mroz %>% filter(!is.na(lwage))

# Wooldridge sugere um teste (baseado em Hausman (1978)) de endogeneidade. É um procedimento simples: (i) basta estimar a forma reduzida do regressor cuja endogeneidade se quer testar, (ii) salvar os resíduos, (iii) estimar o modelo de interesse por OLS incluindo tais resíduos como variável independente e (iv) conduzir um teste de hipótese, preferencialmente robusto, sobre os resíduos. A hipótese nula do teste é de que o regressor estudado é exógeno.
mroz$resíduos_1S = mod_1S[["residuals"]]

modAuxiliar = lm(
  lwage ~ educ + exper + expersq + resíduos_1S,
  data = mroz 
)
coeftest(modAuxiliar, vcov. = vcovHC)

# Observe que este teste, sem robustez, reproduz o "teste de Wu-Hausman", no output pronto do AER (lá se aponta um teste F; basta tomar o quadrado da estatística t daqui). Novamente, a vantagem de fazer a implementação manual é poder conduzir o teste robusto.
summary(mod_2SLS, diagnostics = TRUE)
(coeftest(modAuxiliar)[nrow(coeftest(modAuxiliar)),3])^2

# OBSERVAÇÃO. No caso de múltiplos regressores endógenos, estime uma regressão auxiliar incluindo resíduos para cada um destes regressores, depois, conduza um teste de exclusão conjunta.

# ----------------------------------------------------------------------
# TESTE DE RESTRIÇÕES SOBREIDENTIFICADAS
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('mod_2SLS', 'mod_1S', 'mod_2S')])
data(mroz)
mroz = mroz %>% filter(!is.na(lwage))

# Basta: (i) estimar a equação estrutural por 2SLS, (ii) salvar os respectivos resíduos, (iii) regredir os resíduos contra todas as variáveis exógenas, (iv) salvar o R2 dos resultados e (v) comparar a estatística de teste, como descrita abaixo, com a respectiva distribuição. A hipótese nula do teste é de que todos os instrumentos são, de fato, exógenos.
mroz$resíduos_2SLS = mod_2SLS[['residuals']]

modAuxiliar = lm(
  resíduos_2SLS ~ exper + expersq + motheduc + fatheduc,
  data = mroz
)

R2 = summary(modAuxiliar)[['r.squared']]

estatísticaTeste = nrow(mod_2SLS[['model']]) * R2

# Observe que a estatística de teste tem distribuição qui-quadrado com graus de liberdade iguais ao número de instrumentos (excluídos do modelo) menos o número de variáveis endógenas (no caso, 1).
pchisq(estatísticaTeste, df = 1, lower.tail = FALSE)
