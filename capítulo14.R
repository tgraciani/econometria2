rm(list = ls())
#options(scipen = 100, digits = 7)
source('~/Dropbox/doutorado/monitoria/monitorias/display.R')

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)
library(plm) # Atenção!

### SUMÁRIO ###
# 1 -- ESTIMAÇÃO DE EFEITOS FIXOS (FE)
# 2 -- FE: ESTIMAÇÃO MANUAL
# 3 -- FE: BETWEEN ESTIMATOR
# 4 -- FE: DUMMY VARIABLE REGRESSION
# 5 -- FE: COMPARAÇÃO COM FD
# 6 -- ESTIMAÇÃO DE EFEITOS ALEATÓRIOS (RE)
# 7 -- RE: ESTIMAÇÃO MANUAL
# 8 -- TESTE DE HAUSMAN
# 9 -- ESTIMAÇÃO DE EFEITOS ALEATÓRIOS CORRELACIONADOS (CRE)
# 10 -- CLUSTERING

# Consulte <https://cran.r-project.org/web/packages/wooldridge/wooldridge.pdf> para obter um dicionário das bases de dados do Wooldridge.

# ----------------------------------------------------------------------
# ESTIMAÇÃO DE EFEITOS FIXOS ('PLM')
# ----------------------------------------------------------------------
data(wagepan)

# nr: person identifier
# year: 1980 to 1987
# married: =1 if married 
# educ: years of schooling
# union: =1 if in union
# lwage: log(wage)

# Seleciona as variáveis relevantes.
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)

# Cria dummies para os anos.
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)

# É importante "converter" a base para painel (no formato do 'plm') após fazer modificações como as acima (se você as fizer depois, o 'plm' talvez não consiga estimar o modelo).
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))
pdim(wagepan)

# Estimação de um modelo básico de efeitos fixos via 'plm'. Veja que só se deve incluir as variáveis que aparecem na "within transformation": se você incluir algum regressor sem variabilidade no tempo, o 'plm' o descarta automaticamente.
modFE = plm(
  lwage ~ married + union + (year_1981 + year_1982 + year_1983 + year_1984 + year_1985 + year_1986 + year_1987)*educ,
  #lwage ~ married + union + factor(year)*educ,
  data = wagepan,
  model = 'within'
)
summary(modFE)

# ----------------------------------------------------------------------
# FE: ESTIMAÇÃO MANUAL
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'show_math')])
data(wagepan) # Mesmo conjunto de dados, veja o dicionário acima.

# Seleciona as variáveis relevantes.
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)

# Cria dummies para os anos.
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)

# Cria interações entre dummies de anos e 'educ' (será necessário porque vou fazer as contas matricialmente, ou seja, sem o 'lm'; é preciso então criar todas as variáveis usadas como "colunas independentes").
wagepan = wagepan %>% mutate(
  year_1981educ = year_1981 * educ,
  year_1982educ = year_1982 * educ,
  year_1983educ = year_1983 * educ,
  year_1984educ = year_1984 * educ,
  year_1985educ = year_1985 * educ,
  year_1986educ = year_1986 * educ,
  year_1987educ = year_1987 * educ
)

# Faz o procedimento de "time demeaning". Observe que a operação é feita após se "agrupar" os dados segundo a variável 'nr'. Todas as variáveis no banco de dados (exceto por 'nr') passam pela operação. Note que como "destruímos" as variáveis originais, preferi fazer as mudanças associando-as a um novo objeto.
d = wagepan %>%
  group_by(nr) %>% # Faz o agrupamento adequado.
  mutate_at(
    vars(-group_cols()), # Seleciona as variáveis que passarão pelo demeaning.
    list(~. - mean(.)) # Aplica a operação.
  ) %>%
  ungroup() # Desagrupa os dados.

# Cria as matrizes 'y' e 'X'. Usarei as "fórmulas matriciais" de OLS. Veja que, na definição de 'X', tive que excluir nominalmente todos os regressores que não são regressores no modelo de FE.
y = data.matrix(
  d %>% select(lwage)
)
dim(y)

X = data.matrix(
  d %>% select(-c(nr, lwage, educ, year))
)
dim(X)

# Obtém os coeficientes da estimação por OLS. Observe que a operação '%*%' denota multiplicação de matrizes (como aprendemos em álgebra linear) e que 'solve()' determina a inversa de seu argumento, da forma como usamos abaixo.
b = as.vector(solve(t(X)%*%X) %*% t(X)%*%y)
names(b) = colnames(X)
b

# Comparação com os "estimadores manuais" e os "automatizados".
cbind(
  'Manual' = b,
  'Automático' = coef(modFE),
  'Diferença' = b - coef(modFE)
)

# Calcula os graus de liberdade da estimação (usando o número de indivíduos, de anos e de regressores).
dimN = length(unique(d$nr))
dimT = length(unique(d$year))
dimK = length(b)

# Computa os resíduos.
e = as.vector(y - (X %*% b))

# Calcula a matriz de variância-covariância dos estimadores.
VCV = as.matrix(
  (1 / (dimN * (dimT - 1) - dimK)) * as.numeric((t(e)%*%e)) * solve(t(X)%*%X)
)

# Calcula os erros-padrão (corretos enfim) dos estimadores de OLS do modelo de FE.
SE = as.vector(sqrt(diag(VCV)))

# Exibe nosso output caseiro.
cbind(
  'Estimativa' = b,
  'Erro-Padrão' = SE,
  't' = b/SE,
  'Pr(>|t|)' = 2 * pt(abs(b)/SE, df = (dimN * (dimT - 1) - dimK), lower.tail = FALSE)
)

summary(modFE)

# ----------------------------------------------------------------------
# FE: BETWEEN ESTIMATOR
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'show_math')])
data(wagepan) # Mesmo conjunto de dados, veja o dicionário acima.

# Mesmo procedimento, feio já duas vezes acima.
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))

# Por definição, os "between estimators" são aqueles obtidos via OLS do modelo que se subtrai do populacional para chegar ao "within". Observe que todas as dummies de ano colapsam para um intercepto único.
show_math(
  '$$\\overline{y}_i = \\beta_0 + \\beta_1 \\overline{x}_{1,i} + \\cdots + \\beta_k \\overline{x}_{k,i} + \\overline{\\varepsilon}_{i}$$'
)
modFE_Between = plm(
  lwage ~ married + union + (year_1981 + year_1982 + year_1983 + year_1984 + year_1985 + year_1986 + year_1987)*educ,
  #lwage ~ married + union + educ,
  data = wagepan,
  model = 'between'
)
summary(modFE_Between)

# ----------------------------------------------------------------------
# FE: DUMMY VARIABLE REGRESSION
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'show_math')])
data(wagepan) # Mesmo conjunto de dados, veja o dicionário acima.

# A ideia é estimar o efeito não observado para cada observação na amostra. Implementa-se a estimação pela inclusão de variáveis dummy para cada indivíduo amostrado (na prática, criando N interceptos "particulares").
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)
wagepan = wagepan %>% dummy_cols(select_columns = 'nr', remove_first_dummy = TRUE)
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))

# Para a estimação em si, seria mais simples usar o 'factor()', como acima. Aproveito a seção para apresentar novamente (como fiz em uma lista passada) a construção da fórmula do 'lm' me valendo dos padrões dos nomes das variáveis e da função 'paste()'. (Se você decidir fazer com o 'factor()', precisará fazer ajustes nas tabelas de comparação mais abaixo.)
show_math(
  '$$\\mathbf{y} = \\mathbf{X} \\beta + \\mathbf{d}_i \\lambda + \\varepsilon$$'
)
modDVR = lm(
  as.formula(
    paste(
      'lwage ~ married + union',
      paste(grep('nr_', names(wagepan), value = TRUE), collapse = ' + '),
      paste(grep('year_', names(wagepan), value = TRUE), collapse = ' + '),
      paste(grep('year_', names(wagepan), value = TRUE), 'educ', sep = ':', collapse = ' + '),
      sep = ' + '
    )
  ),
  data = wagepan
)
summary(modDVR)

# Salva coeficientes e erros-padrão da estimação.
bDVR = coef(modDVR)
seDVR = coef(summary(modDVR))[, 'Std. Error']

# Retoma as mesmas informações, mas do modelo de FE.
bFE = coef(modFE)
seFE = coef(summary(modFE))[, 'Std. Error']

# Verificando a previsão de que os estimadores de FE e DV são idênticos (logicamente, para os regressores presentes no modelo de FE).
cbind(
  'FE' = bFE,
  'DVR' = bDVR[names(bFE)]
)

# Erros-padrão também são iguais.
cbind(
  'FE' = seFE,
  'DVR' = seDVR[names(bFE)]
)

# É possível, ainda, testar conjuntamente a presença de heterogeneidades individuais (mas veja que isso não significa que os componentes não observáveis são correlacionados com os regressores do modelo).
# OBSERVAÇÃO: O comando abaixo está correto, mas não consegui executá-lo (ele bate em algum limite do meu computador que não consegui flexibilizar). Talvez funcione com quem tiver um computador mais poderoso...
#seleçãoTeste = grep('nr_', names(wagepan), value = TRUE)
#linearHypothesis(modDVR, seleçãoTeste)

# Uma alternativa é calcular a estatística F "na mão". Para isso, precisamos da soma dos quadrados dos resíduos dos modelos, irrestrito e restrito, e os respectivos graus de liberdade (vide capítulo 4 do Wooldridge).
show_math(
  '$$\\text{F} = \\frac{(\\text{SSR}_r-\\text{SSR}_{ur})/q}{\\text{SSR}_{ur}/(n-k-1)}$$'
)

SSR_ur = sum(modDVR[['residuals']]^2)

modDVR_r = lm(
  as.formula(
    paste(
      'lwage ~ married + union',
      paste(grep('year_', names(wagepan), value = TRUE), collapse = ' + '),
      paste(grep('year_', names(wagepan), value = TRUE), 'educ', sep = ':', collapse = ' + '),
      sep = ' + '
    )
  ),
  data = wagepan
)
SSR_r = sum(modDVR_r[['residuals']]^2)

dfr_num = modDVR_r[['df.residual']] - modDVR[['df.residual']] # length(unique(wagepan$nr))-1
dfr_den = modDVR[['df.residual']]

# Então basta usar a fórmula do teste como vocês viram no primeiro curso de econometria de vocês.
FF = ((SSR_r - SSR_ur) / dfr_num) / (SSR_ur / dfr_den)
FF
pf(FF, dfr_num, dfr_den, lower.tail = FALSE)

# ----------------------------------------------------------------------
# FE: COMPARAÇÃO COM FD
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'show_math')])
data(wagepan)

# Procedimentos similares aos feitos com a mesma base de dados. Esta é restrita a dois anos de dados para obter a equivalência entre FD e FE.  Também se cria uma dummy indicando as informações do segundo período.
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)
wagepan = wagepan %>% filter(year %in% 1980:1981)
wagepan = wagepan %>% mutate(year_81 = if_else(year == 1980, 0, 1, missing = NA_real_))
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))

# Estimadores de FE e FD coincidem quando (i) só se têm dois períodos de dados e (ii) inclui-se uma dummy para o segundo período na especificação de FE (cuja estimativa bate com o intercepto de FD).
modFE_2 = plm(
  lwage ~ year_81 + married + union + year_81:educ,
  data = wagepan,
  model = 'within'
)
summary(modFE_2)

modFD = plm(
  lwage ~ married + union + year_81:educ,
  # lwage ~ married + union + year_81 + year_81:educ,
  data = wagepan,
  model = 'fd'
)
summary(modFD)

# ----------------------------------------------------------------------
# ESTIMAÇÃO DE EFEITOS ALEATÓRIOS (RE)
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'show_math')])
data(wagepan)

show_math(
  '$$y_{i,t} = a_i + \\mathbf{x}_{i,t}^T \\beta + \\varepsilon_{i,t}$$',
  '$$y_{i,t} = \\mathbf{x}_{i,t}^T \\beta + \\nu_{i,t}$$',
  '$$\\text{corr}(\\nu_{i,t},\\nu_{i,s}) = \\frac{\\sigma^2_a}{\\sigma^2_a + \\sigma^2_u}$$'
)

# Procedimento usual.
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))

# A estimação de um modelo de RE é, como nos casos de FD e FE, direta com o pacote 'plm'. Veja que há uma série de opções sobre os estimadores das variâncias necessárias para construir o modelo de RE (detalhes técnicos demais; só notem que todos os estimadores devem ser consistentes, embora possam duvergir sob amostras finitas.)
show_math(
  '$$(y_{i,t} - \\theta \\overline{y}_i) = (\\mathbf{x}_{i,t}^T - \\theta  \\overline{\\mathbf{x}}_i^T) \\beta + (\\nu_{i,t}  - \\theta \\overline{\\nu}_i)$$',
  '$$\\theta = 1 - \\sqrt{\\frac{\\sigma^2_\\varepsilon}{\\sigma^2_\\varepsilon + T \\sigma^2_a}}$$'
)
modRE = plm(
  lwage ~ educ + black + hisp + exper + married + union + year_1981 + year_1982 + year_1983 + year_1984 + year_1985 + year_1986 + year_1987,
  data = wagepan,
  #random.method = 'swar', # Default.
  #random.method = 'amemiya',
  #random.method = 'walthus',
  #random.method = 'nerlove',
  model = 'random'
)
summary(modRE)

# Para efeitos de comparação, é possível extrair o estimador de θ do objeto gerado pelo 'plm'.
modRE[['ercomp']]

# ----------------------------------------------------------------------
# RE: ESTIMAÇÃO MANUAL
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'modRE', 'show_math')])
data(wagepan)

# Procedimentos usuais.
wagepan = wagepan %>% select(nr, lwage, educ, black, hisp, exper, married, union, year)
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)

# Há várias formas de se estimar um modelo de RE manulamnete. Vou seguir o Wooldridge (“Econometric Analysis of Cross Section and Panel Data”, 2010, páginas 291-296). Primeiro, estima-se um modelo linear, via OLS, fazendo pooling dos dados e usando as mesmas variáveis do modelo de RE.
modPOLS = lm(
  lwage ~ educ + black + hisp + exper + married + union + year_1981 + year_1982 + year_1983 + year_1984 + year_1985 + year_1986 + year_1987,
  data = wagepan
)

# Salva uma série de dimensões necessárias nos cálculos.
dimN = length(unique(wagepan$nr))
dimT = length(unique(wagepan$year))
dimK = length(coef(modPOLS))

# Para implementar as etimativas, constrói uma matriz 'V' com os resíduos de POLS. Em particular, linhas indicam diferentes anos de amostra e colunas, indivíduos.
V = matrix(modPOLS[['residuals']], nrow = dimT)
dim(V)

# Na notação do Wooldrige, é necessiário estimar "σ^2_e" (variância do modelo populacional) e "σ^2_a" (variância do efeito não observado). É mais fácil, contudo, (i) estimar a variância do erro composto (σ^2_v), (ii) estimar "σ^2_a" diretamente e (iii) obter "σ^2_e" como um resíduo.
show_math(
  '$$\\sigma^2_{\\nu} = \\sigma^2_{\\varepsilon} + \\sigma^2_a$$'
)

#A parte (i) é a mais simples:
sigma2v = 1/(dimN * dimT - dimK) * sum(modPOLS[['residuals']]^2)

# Forma alternativa de se implementar (i), similar àquela que vou apresentar para (ii):
show_math(
  '$$\\hat{\\sigma}^2_\\nu = \\frac{1}{NT-K} \\sum_{i=1}^{N} \\sum_{t=1}^{T} \\hat{\\varepsilon}_{i,t}^2 $$'
)
SOMA = 0
for (i in 1:ncol(V)) {
  for(t in 1:nrow(V)) {
    SOMA = SOMA + (V[t,i]^2)
  }
}
sigma2v = 1/(dimN * dimT - dimK) * SOMA

# Wooldridge sugeriu estimar "σ^2_a" usando três somatórios encadeados. Como a estrutura da conta não é tão simples, achei mais fácil (embora certamente seja computacionalmente ineficiente) fazer a estimação dentro de um loop. Em particular, aqui estou seguindo a equação 10.37, na página 296.
show_math(
  '$$\\hat{\\sigma}^2_a = \\frac{1}{NT(T-1)2^{-1}-K} \\sum_{i=1}^{N} \\sum_{t=1}^{T-1} \\sum_{s=t+1}^{T} \\hat{\\varepsilon}_{i,t} \\hat{\\varepsilon}_{i,s} $$',
  css = 'color: back; font-size: 20px;'
)
SOMA = 0
for (i in 1:ncol(V)) {
  for (t in 1:(nrow(V)-1)) {
    for (s in (t+1):nrow(V)) {
      SOMA = SOMA + (V[t,i] * V[s,i])
    }
  }
}
sigma2a = (1/((dimN*dimT*(dimT-1)/2 - dimK))) * SOMA

# Por fim, "σ^2_e" sai como um resíduo:
sigma2e = sigma2v - sigma2a

# Uma vez munidos de "σ^2_a" e "σ^2_e", para estimar o θ basta seguir o Wooldridge (do curso mesmo), como descrito na página 470 (estou consultando a sétima edição).
theta = 1 - sqrt(sigma2e / (sigma2e + dimT * sigma2a))
theta

# Para efeitos de comparação, observe que há uma diferença (ainda que pequena) entre o cálculo manual sugerido no livro-texto e o implementado pelo 'plm'. Possivelmente esta decorre dos estimadores escolhidos para as variâncias necessárias.
theta
modRE[['ercomp']]['theta']

# Com uma estimativa de θ, resta fazer o demeaning dos dados. Observe, contudo, que é preciso fazer o procedimento também para o intercepto. Por isso crio uma variável 'Intercept' igual a um (como já mencionei algumas vezes, inserir um intercepto na estimação significa considerar a inclusão de um regressor igual a 1 para todas as observações).
wagepan$Intercept = 1

# Imlementa o demeaning de RE. Observe que o procedimento é feito sobre todas as variáveis (exceto 'nr'), então é recomendável fazer o assigment para um objeto novo.
d = wagepan %>%
  group_by(nr) %>%
  mutate_at(
    vars(-group_cols()),
    list(~. - theta*mean(.))
  ) %>%
  ungroup()

# Por fim, obtêm-se os estimadores de RE estimando o modelo por OLS com as variáveis transformadas. Note que o '0' na fórmula indica para o 'lm' que queremos uma regressão sem intercepto (porque estamos provendo o nosso próprio, transformado).
modRE_2 = lm(
  lwage ~ 0 + Intercept + educ + black + hisp + exper + married + union + year_1981 + year_1982 + year_1983 + year_1984 + year_1985 + year_1986 + year_1987,
  data = d
)
summary(modRE_2)

# Tabela para comparação. Veja que a diferença segue das diferentes estimativas de θ.
cbind(
  'PLM/Swar' = coef(modRE),
  'Manual' = coef(modRE_2),
  'Dif.' = coef(modRE) - coef(modRE_2)
)

# ----------------------------------------------------------------------
# TESTE DE HAUSMAN
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'modRE', 'show_math')])

# Há um comando pronto para executar o teste de Hausman para comparar modelos de FE e RE. De modo geral, o teste verifica se dois estimadores são arbitrariamente próximos (iguais). No caso de a hipótese nula ser rejeitada, é costumeiro ficar com os de FE (embora, estatisticamente falando, o teste não aponte nenhuma vantagem de FE sobre RE).
# OBSERVAÇÃO. Apesar de o "teste rodar" sem avisos, Wooldridge (2010, p. 329) discute a não admissibilidade da inclusão de (i) variáveis constantes no tempo (naturalmente excluídas do modelo de FE, mas não do RE) e (ii) variáveis de "aggregate time effects", que só variam na dimensão do tempo (como, por exemplo, dummies de tempo; o problema é "avançado" para o que vocês estão vendo, este é só um comentário "extra"). É importante estar atento que o comando do 'plm' só vai comparar no teste coeficientes comuns a ambos os modelos.
# A distribuição do teste é derivada sob a hipótese de que RE é consistente (ie, que as hipóteses básicas do método, listadas no Wooldridge, são válidas).
phtest(modFE, modRE)

# ----------------------------------------------------------------------
# ESTIMAÇÃO DE EFEITOS ALEATÓRIOS CORRELACIONADOS (CRE)
# ----------------------------------------------------------------------
rm(list = ls()[!ls() %in% c('modFE', 'modRE', 'show_math')])
data(wagepan) # Mesmo conjunto de dados, veja o dicionário acima.

# Prodedimentos usuais.
wagepan = wagepan %>% select(nr, lwage, married, union, year, educ)
wagepan = wagepan %>% dummy_cols(select_columns = 'year', remove_first_dummy = TRUE)
wagepan = wagepan %>% mutate(
  year_1981educ = year_1981 * educ,
  year_1982educ = year_1982 * educ,
  year_1983educ = year_1983 * educ,
  year_1984educ = year_1984 * educ,
  year_1985educ = year_1985 * educ,
  year_1986educ = year_1986 * educ,
  year_1987educ = year_1987 * educ
)
wagepan = pdata.frame(wagepan, index = c('nr', 'year'))

# Também há no 'plm' uma implementação automatizada do do modelo CRE uma vez que permite a estimação de RE e disponibiliza a função 'Between()' que, lembre-se, calcula médias temporais (restritas a cada indivíduo).
# Observe que, como o painel usado é balanceado, não é preciso adicionar a função 'Between()' para as dummies de tempo (desde que se permita um intercepto) nem para a interação de outras variáveis com estas dummies (desde que se inclua a variável em si, geralmente que não se acomodam no modelo de FE.
modCRE = plm(
  lwage ~ married + union + year*educ + Between(married) + Between(union) + Between(year*educ),
  data = wagepan,
  model = 'random'
)

# Observe que só se garante a coincidência dos estimadores de FE e CRE se se incluir o 'Between()' de todos os regressores presentes no modelo de FE (o que não foi feito na equação) 
summary(modCRE)
summary(modFE)

# Esta especificação permite a condução de um teste para a escolha entre os modelos de FE e RE: trata-se de um simples teste F sobre as variáveis adicionadas via 'Between()'. Sua hipótese nula é de que não se precisa do modelo CRE ou FE, ou seja, que RE são suficientes.
linearHypothesis(modCRE, matchCoefs(modCRE, 'Between'))

# ----------------------------------------------------------------------
# CLUSTERING
# ----------------------------------------------------------------------
# As observações do script anterior sobre clustering continuam válidas aqui: para quaisquer modelos estimados, você pode usar opções de clustering com os comandos 'coeftest' ou 'linearHypothesis' do mesmo modo como feito antes.
