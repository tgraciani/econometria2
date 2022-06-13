rm(list = ls())

library(wooldridge)
library(tidyverse)
library(fastDummies)
library(lmtest)
library(car)
library(plm) # Atenção!

### SUMÁRIO ###
# 1 -- EMPILHAMENTO DE BASES: 'RBIND'
# 2 -- "MERGE" DE BASES: 'CBIND' E 'MERGE'
# 3 -- MODELOS LINEARES E POOLED CROSS-SECTIONS (POLS)
# 4 -- TESTE DE CHOW
# 5 -- DIFFERENCES-IN-DIFFERENCES ESTIMATOR (DD)
# 6 -- CONSTRUÇÃO DE PAINÉIS: PACOTE 'PLM'
# 7 -- OPERAÇÕES COM PAINÉIS
# 8 -- FIRST-DIFFERENCE ESTIMATOR (FD)
# 9 -- FD: AUTOCORRELAÇÃO SERIAL E CLUSTERING

# Consulte <https://cran.r-project.org/web/packages/wooldridge/wooldridge.pdf> para obter um dicionário das bases de dados do Wooldridge.

# ----------------------------------------------------------------------
# EMPILHAMENTO DE BASES: 'RBIND'
# ----------------------------------------------------------------------
# Cria dois bancos de dados hipotéticos (usando extrações aleatórias de diferentes distribuições). Observe como a ordem das variáveis é diferente entre os objetos (mas isso não é um problema).
d1 = as.data.frame(
  matrix(
    c(
      1:10,
      rnorm(10, mean = 7, sd = 3),
      rexp(10, rate = 5),
      runif(10, min = -10, max = 10)
    ),
    nrow = 10,
    ncol = 4
    )
  )
colnames(d1) = c('ID', 'y', 'x1', 'x2')

d2 = as.data.frame(
  matrix(
    c(
      runif(10, min = -10, max = 10),
      rnorm(10, mean = 7, sd = 3),
      rexp(10, rate = 5),
      101:110
    ),
    nrow = 10,
    ncol = 4
    )
  )
colnames(d2) = c('x2', 'y', 'x1', 'ID')

# Para empilhar as bases de dados ("encaixando" uma embaixo da outra), use o comando 'rbind'.
d = rbind(d1, d2)

# Eventualmente pode ser necessário ou útil rearranjar as colunas de um banco de dados de modo compatível com a ordem de outro. No exemplo abaixo, reorganizam-se as colunas de 'd2' seguindo a ordem de 'd1'.
d3 = d2[,names(d1)]
rm(d3)

# Como observação final, observe que, se não houver coincidência das colunas de ambas as bases, o 'rbind' retorna um erro. Nestes casos, você deve trabalhar para compatibilizar as bases.
dd2 = d2[c('x2', 'y', 'ID')]
d = rbind(d1, dd2)

# A primeira solução consiste em só considerar as colunas da "menor base", com menos variáveis.
d = rbind(d1[names(dd2)], dd2)

# A segunda é preencher a base menor com colunas de NA. (Apesar de funcionar, veja que a primeira solução é computacionalmente superior: não é de serventia alguma manter variáveis que não poderão ser usadas nas estimações).
dd2$x1 = NA
d = rbind(d1, dd2)
rm(dd2)

# ----------------------------------------------------------------------
# "MERGE" DE BASES: 'CBIND' E 'MERGE'
# ----------------------------------------------------------------------
# Muda a variável de 'ID' (supondo que temos duas entrevistas com as mesmas pessoas em dois momentos diferentes).
d2$ID = 1:10

# Renomeia as variáveis do bloco anterior. Interprete 'xij' como a variável 'i' no tempo 'j'.
colnames(d1) = row.names = c('ID', 'y1', 'x11', 'x21')
colnames(d2) = row.names = c('x22', 'y2', 'x12', 'ID')

# Se, ao invés de "empilhar" as bases você precisar justapô-las lado a lado, use o 'cbind' (como tenho feito com as "tabelinhas" nos últimos scripts).
d = cbind(d1, d2)
#d = cbind(d1, d2[,!(names(d2) %in% c('ID'))])

# Entretanto, o 'cbind' pode ser problemático se as observações não estiverem ordenadas (nos exemplos aqui, a variável 'ID'). No script abaixo, se aleatoriza a ordem das linhas de 'd2' e se faz um 'cbind' problemático.
dd2 = d2[sample(nrow(d2)),]
d = cbind(d1, dd2)
rm(dd2)

# A solução, nestes casos, é o 'merge'. Basta especificar qual a variável identifica as observações em cada base de dados.
d = merge(d1, d2, by = 'ID')

# Analogamente ao empilhamento, diferença nos tamanhos das linhas gera problemas com o 'cbind', mas não com o 'merge' (que descarta observações não presentes em ambas as bases).
dd2 = d2[c(1:6, 8:10),]
cbind(d1, dd2)
d = merge(d1, dd2, by = 'ID')

# ----------------------------------------------------------------------
# MODELOS LINEARES E POOLED CROSS-SECTIONS (POLS)
# ----------------------------------------------------------------------
rm(list = ls())
data('cps78_85')

# educ: years of schooling
# female: =1 if female
# exper: age - educ - 6
# union: =1 if belong to union
# lwage: log hourly wage
# y85: =1 if year == 85

# A implementação de POLS é idêntica à de OLS, só varia mesmo a "natureza" da base de dados. Só observe que como se usa o logaritmo dos salários como variável dependente, a inclusão de um intercepto para o segundo período nos livra da necessidade de corrigir a variável pela inflação do período.
# Veja que na sintaxe da fórmula do 'lm', 'x1 * x2' significa 'x1 + x2 + x1:x2', onde este último termo é a interação.
modPOLS = lm(
  lwage ~ y85 * (educ + female) + exper + I((exper^2) / 100) + union,
  #lwage ~ y85 + educ + female + y85:educ + y85:female + exper + I((exper^2) / 100) + union,
  data = cps78_85)
summary(modPOLS)

# Um ponto interessante é o de que a inclusão de todas as possíveis interações com a dummy de tempo é equivalente à estimação dois modelos separadamente (relativamente aos coeficientes), cada um com dados de um período.
mod78 = lm(
  lwage ~ educ + female + exper + I((exper^2) / 100) + union,
  data = cps78_85,
  subset = c(year == 78))

mod85 = lm(
  lwage ~ educ + female + exper + I((exper^2) / 100) + union,
  data = cps78_85,
  subset = c(year == 85))

mod78_85 = lm(
  lwage ~ y85 * (educ + female + exper + I((exper^2) / 100) + union),
  data = cps78_85)

b = coef(mod78_85)
bBase = b[c(1, 3:7)]
bInteração = b[c(2, 8:12)]

cbind(
  'OLS 1978' = mod78[['coefficients']],
  'OLS 1985' = mod85[['coefficients']],
  'Base {A}' = bBase,
  'Inter. {B}' = bInteração,
  '{A} + {B}' = bBase + bInteração
)

# ----------------------------------------------------------------------
# TESTE DE CHOW
# ----------------------------------------------------------------------
rm(list = setdiff(ls(), 'mod78_85'))

# Para implementar um teste que permite correção para problemas de heterocedasticidade, é necessário usar uma regressão "completa", com dados de todos os períodos e interações com dummies de ano. No caso, optei por manter a mudança de intercepto fora do teste (é um ponto contextual).
# Detalhes (e referências) das opções de 'hccm' estão disponíveis na documentação de 'car'.
linearHypothesis(
  mod78_85,
  paste('y85', c('educ', 'female', 'exper', 'I((exper^2)/100)', 'union'), sep = ':'),
  #c(y85:educ, y85:female, y85:exper, y85:I((exper^2)/100), y85:union),
  vcov. = hccm
  #vcov. = hccm(mod78_85, type = 'hc0')
  #vcov. = hccm(mod78_85, type = 'hc1')
  #vcov. = hccm(mod78_85, type = 'hc2')
  #vcov. = hccm(mod78_85, type = 'hc3') # Default.
  #vcov. = hccm(mod78_85, type = 'hc4')
)

# ----------------------------------------------------------------------
# DIFFERENCES-IN-DIFFERENCES ESTIMATOR (DD)
# ----------------------------------------------------------------------
rm(list = ls())
data('kielmc')

# age: age of house
# rooms: # rooms in house
# area: square footage of house
# land: square footage lot
# baths: # bathrooms
# y81: =1 if year == 1981
# nearinc: =1 if dist <= 15840
# rprice: price, 1978 dollars

# Quando se possoi dados experimentais ou quase-experimentais, a abordagem "padrão" é a de DD. Se não houver regressores adicionais (além da dummy para tratamento), uma forma de obter o efeito do tratamento consiste na estimação de dois modelos, para os períodos pré e pós-tratamento. O estimador de DD é a diferença entre os coeficientes dos dois modelos (associado àquela dummy).
mod78 = lm(
  rprice ~ nearinc,
  data = kielmc,
  subset = (year == 1978))
summary(mod78)
bPré = as.numeric(coef(mod78)[2])
bPré

mod81 = lm(
  rprice ~ nearinc,
  data = kielmc,
  subset = (year == 1981))
summary(mod81)
bPós = as.numeric(coef(mod81)[2])
bPós

DD = bPós - bPré
DD

# Uma desvantagem imediata é não ter um cálculo direto para o erro-padrão do estimador. Para calcular estes mais facilmente (permitindo, inclusive, a correção de qualquer heterocedasticidade), pode-se recorrer a uma única estimação via OLS. No caso, interessa-se pelo coeficiente de interação entre tratamento e dummy de tempo.
modDD = lm(
  rprice ~ nearinc * y81,
  # rprice ~ nearinc + y81 + nearinc:y81
  data = kielmc)
summary(modDD)

# Verifique que as estimativas batem.
all.equal(
  as.numeric(modDD[['coefficients']][4]), DD
)

# Faz teste (robusto) sobre o coeficiente de interesse.
coeftest(mod, vcov. = hccm)

# Como observação final, a inclusão de controles finda a possibilidade de obter a estimativa de DD usando duas regressões separadamente. Nestes casos é necessário obter o estimador de interesse via uma única regressão.
modDD = lm(
  rprice ~ (nearinc * y81) + age + I(age^2) + log(intst) + log(land) + log(area) + rooms + baths,
  data = kielmc)
summary(mod)
bDD = as.numeric(coef(mod)[11])

mod78 = lm(
  rprice ~ nearinc + age + I(age^2) + log(intst) + log(land) + log(area) + rooms + baths,
  data = kielmc,
  subset = (year == 1978))
bPré = as.numeric(coef(mod78)[2])

mod81 = lm(
  rprice ~ nearinc + age + I(age^2) + log(intst) + log(land) + log(area) + rooms + baths,
  data = kielmc,
  subset = (year == 1981))
bPós = as.numeric(coef(mod81)[2])

# Calcula (incorretamente) a estimativa de DD.
bPós - bPré

# Comparação com a correta, 'b[11]'.
all.equal(bDD, bPós - bPré)

# Como observação final, note que a estimação de DD admite a "parallel trends assumption" (qualquer tendência na variável dependente é igual nos grupos de tratados e controles, ou seja, o tratamento não a afeta).
# Uma flexibilização para "contornar" a hipótese consiste no uso do estimador de tripla diferença (DDD). Este, contudo, requer informações de um segundo grupo de controle. (Veja detalhes no Wooldridge, não vou fazer a implementação: é só mais um modelo linear estimado por OLS, só sendo ajustada a especificação.)

# ----------------------------------------------------------------------
# CONSTRUÇÃO DE PAINÉIS: PACOTE 'PLM'
# ----------------------------------------------------------------------
rm(list = ls())
data('crime2')

# pop: population
# crimes: total number index crimes
# unem: unemployment rate
# year: 82 or 87
# ccrmrte: change in crmrte

# É possível trabalhar com painéis sem a ajuda do 'plm', mas é necessária muita cautela para não cometer "gafes" (e.g., diferenciar uma variável no tempo para dois indivíduos diferentes). Em princípio, trabalhem com o 'lpm'. Uma vez confortáveis, experimentem fazer o trabalho "na mão".
# Veja que os dados deste exemplo vêm formatados de um modo específico (cada região, nosso "i", tem duas linhas apresentadas sequencialmente, uma com informações de 1982 e outra de 1985, o "t").
crime2[1:5, c('year', 'pop', 'crimes', 'crmrte', 'unem')]

# Veja o arquivo de ajuda sobre a função 'pdata.frame' para entender as possibilidades acomodadas pela opção 'index'. Abaixo vou explorar outros contextos. Neste caso agora, indica-se o número de indivíduos (requer que os dados estejam bem organizados, como indicado acima).
crime2 = pdata.frame(crime2, index = 46)
crime2[1:5, c('id', 'time', 'year', 'pop', 'crimes', 'crmrte', 'unem')]

# Lista as dimensões da base.
pdim(crime2)

# Uma vez convertida a base para o "tipo" 'pdata.frame', é possível usar os recursos do pacote 'lpm', como operações (para criação de novas variáveis) e estimação dos "modelos canônicos" de dados em painel.

# ----------------------------------------------------------------------
# OPERAÇÕES COM PAINÉIS
# ----------------------------------------------------------------------
rm(list = ls())
data(crime4)

# county: county identifier
# year: 81 to 87
# crmrte: crimes committed per person

# Transforma os dados num painel compatível com o 'plm'. Aqui, diferentemente do outro exemplo, se especificam as variáveis que identificam indivíduos ("i") e tempo ("t"). Não requer tanta organização como antes.
crime4 = pdata.frame(crime4, index = c('county', 'year'))

# Faz uma série de operações sobre a variável 'crmrte': (1) calcula sua defasagem, (2) toma a primeira diferença, (3) computa a média dela considerando cada observação individualmente (agrupando os dados nesta dimensão) e (4) deduz a variável para cada período como desvios em relação à média (novamente, agrupando por observação).
crime4$op1 = lag(crime4$crmrte)
crime4$op2 = diff(crime4$crmrte)
crime4$op3 = Between(crime4$crmrte)
crime4$op4 = Within(crime4$crmrte)
crime4[1:10, c('county', 'year', 'crmrte', 'op1', 'op2', 'op3', 'op4')]

# Infelizmente, o 'dplyr' (pacote chamado pelo 'tidyverse') não funciona com o 'plm'. Existe uma opção, o 'pmdplyr', para quem quiser se aventurar (não vou apresentar, não o conheço).

# ----------------------------------------------------------------------
# FIRST-DIFFERENCE ESTIMATOR (FD)
# ----------------------------------------------------------------------
rm(list = ls())
data(crime2)

# unem: unemployment rate
# crmrte: crimes per 1000 people

# Converte a base de dados para um formato compatível com o 'plm' (é a mesma de outro exemplo acima).
crime2 = pdata.frame(crime2, index = 46)

# O pacote 'plm' possui um comando análogo ao 'lm', permitindo a estimação de todos os modelos adequados a dados em painel deste curso. Para especificar o estimador de primeira diferença, basta usar a opção 'fd'. Observe que, no output, o intercepto se refere a dummy de tempo no "modelo original", antes de se tomar a primeira diferença. (Este ponto do intercepto fica mais claro na estimação manual, abaixo.)
modFD = plm(
  crmrte ~ unem,
  # crmrte ~ d87 + unem,
  data = crime2,
  model = 'fd'
)
summary(modFD)

# É simples implementar manualmente o cômputo da estimativa de primeira diferença sem o comando específico do 'plm'. Só é necessário criar novas variáveis com as primeiras diferenças e estimar o modelo linear por OLS. (Mas, como feito abaixo, é necessário usar o 'plm' para organizar os dados.)
crime2$dCrmrte = diff(crime2$crmrte)
crime2$dUnem = diff(crime2$unem)
crime2$dd87 = diff(crime2$d87)

modFD = lm(
  dCrmrte ~ dUnem,
  #dCrmrte ~ 0 + dd87 + dUnem,
  data = crime2
)
summary(modFD)

# ----------------------------------------------------------------------
# FD: AUTOCORRELAÇÃO SERIAL E CLUSTERING
# ----------------------------------------------------------------------
rm(list = ls())
data(ezunem)
ezunem = pdata.frame(ezunem, index = 22)

# year: 1980 to 1988
# uclms: unemployment claims
# ez: =1 if have enterprise zone
# luclms: log(uclms)

# Obtém estimadores de FD para discussão abaixo.
modFD = plm(
  luclms ~ ez + d81 + d82 + d83 + d84 + d85 + d86 + d87 + d88,
  data = ezunem,
  model = 'fd'
)
summary(modFD)

# O Wooldridge sugere um teste de autocorrelação serial baseado nos resíduos da estimação (é um modelo simples de séries temporais).
u = modFD[['residuals']]
u = u[match(rownames(ezunem), names(u))]
ezunem$u = u
ezunem$Lu = lag(ezunem$u)

modTeste = lm(
  u ~ 0 + Lu,
  data = ezunem
)
summary(modTeste)

# Havendo evidência de autocorrelação no teste acima, é recomendado calcular "cluster-robust standard errors". Nós não vamos discutir os detalhes técnicos (são mais complicados que os da estimação robusta que já tinhamos usado). Basta estarem cientes de que (i) procedimentos "comuns" de robustez não corrigem autocorrelação serial, (ii) que é um problema "natural" de dados em painel e (iii) cuja solução moderna é o uso de clustering.
# Há vários algoritmos de correção (como com as estimações robustas vistas até agora). Uma opção disponível é a 'sss' ("Stata Small Sample"), calibrada para amostras pequenas.
# Consulte a documentação do pacote 'plm' para obter as referências de cada procedimento. Veja que além de escolher um "tipo" (HC0, sss, HC1, HC2, HC3, HC4), você também pode escolher um "método". No caso, a opção 'arellano' é a única que corrige autocorrelação serial.
coeftest(modFD)
coeftest(modFD, vcovHC(modFD, type = 'sss'))

# A correção também funciona com o 'linearHypothesis' (que vocês podem usar para fazer testes conjuntos).
linearHypothesis(modFD, 'ez', vcov. = vcovHC(modFD, type = 'sss'))
