# Magma_Trace: Simulador de Evolução Geoquímica de Elementos-Traço

[![License: MIT](https://img.shields.io/badge/License-MIT-darkred.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**Autor:** Felipe Padilha Leitzke  
**Afiliação:** Departamento de Mineralogia e Petrologia (DEMIPE), Instituto de Geociências (IGEO), Universidade Federal do Rio Grande do Sul (UFRGS)  
**Ano:** 2026 | **Licença:** MIT

---

## Escopo Científico

O **Magma_Trace** é um aplicativo desktop interativo desenvolvido em Python (`PyQt5`, `NumPy`, `Pandas`, `Matplotlib`) para a modelagem quantitativa do fracionamento de elementos-traço durante processos petrogenéticos ígneos. Voltado tanto ao ensino de graduação e pós-graduação em Geoquímica e Petrologia quanto à exploração paramétrica rápida de modelos magmáticos, o programa resolve analiticamente a evolução de 25 elementos-traço (LILE, HFSE, ETR de La a Lu e metais de transição compatíveis) em sistemas **Basálticos (Máficos)**, **Andesíticos (Intermediários)** e **Riolíticos (Félsicos)**.

---

## Modelos Termodinâmicos e Equacionamento

Nos modelos de anatexia e fusão parcial (Shaw, 1970), o coeficiente de partição global inicial da fonte ($D_0$) e o parâmetro ponderado pela reação de fusão ($P$) são definidos para cada elemento como:

$$D_0 = \sum_{i} X_{0,i} \, K_{d,i} \qquad \text{e} \qquad P = \sum_{i} p_i \, K_{d,i}$$

onde $X_{0,i}$ é a fração modal em peso do mineral $i$ na fonte ($\sum X_{0,i} = 1$) e $p_i$ é a fração estequiométrica líquida com que o mineral $i$ entra no fundido ($\sum p_i = 1$).

* **Fusão Modal ($p_i = X_{0,i}$):** Implica $P = D_0$ e consumo proporcional de todas as fases até $F = 1$.
* **Fusão Não-Modal Cotética ($p_i \ge 0, \, p_i \neq X_{0,i}$):** Todas as fases reagentes são consumidas em taxas distintas de suas modas iniciais.
* **Fusão Não-Modal Peritética / Incongruente ($p_i < 0$):** Permite modelar reações em que uma fase mineral é precipitada enquanto as demais fundem (ex.: fusão de lherzolito a espinélio gerando olivina peritética com $p_{\text{Ol}} < 0$, ou fusão por desidratação de biotita/anfibólio gerando granada ou ortopiroxênio peritético).
* **Limite Físico de Exaustão Mineral ($F_{\max}$):** A massa residual de cada fase no resíduo sólido é $m_i(F) = X_{0,i} - p_i F$. O programa monitora e trunca automaticamente a reação no ponto exato em que a primeira fase consumida ($p_i > 0$) se esgota do resíduo:

$$F_{\max} = \min_{p_i > X_{0,i}} \left( \frac{X_{0,i}}{p_i} \right)$$

### 1. Fusão em Equilíbrio Modal e Não-Modal (*Batch Melting* — Shaw, 1970)
$$\frac{C_l}{C_0} = \frac{1}{D_0 + F(1 - P)}$$

### 2. Fusão Fracionada Modal e Não-Modal (*Fractional Melting* — Shaw, 1970)
* **Líquido Acumulado (*Pooled Fractional Melt*):**
$$\frac{\bar{C}_l}{C_0} = \frac{1}{F} \left[ 1 - \left( 1 - \frac{P \cdot F}{D_0} \right)^{1/P} \right]$$
* **Líquido Instantâneo (*Instantaneous Melt*):** $\frac{C_l}{C_0} = \frac{1}{D_0} \left( 1 - \frac{P \cdot F}{D_0} \right)^{\left(\frac{1}{P} - 1\right)}$ *(Com tratamento analítico exato no limite* $P \to 0$*, onde o termo potencial converge para* $\exp(-F/D_0)$*, e avaliação logarítmica estabilizada via* `expm1`/`log1p`*).*

### 3. Cristalização Fracionada (Destilação de Rayleigh)
$$\frac{C_l}{C_0} = F^{(D - 1)}$$

### 4. Assimilação e Cristalização Fracionada (AFC — DePaolo, 1981)
Sendo $r = \dot{M}_a / \dot{M}_c$ (com $r < 1$) a razão entre a taxa de assimilação crustal e a taxa de cristalização, e $z = \frac{r + D - 1}{r - 1}$:

$$\frac{C_l}{C_0} = F^{-z} + \left( \frac{r}{r + D - 1} \right) \frac{C_a}{C_0} \left( 1 - F^{-z} \right)$$

* **Tratamento analítico da singularidade $z \to 0$ ($r + D = 1$):**
$$\frac{C_l}{C_0} = 1 + \left( \frac{r}{r - 1} \right) \frac{C_a}{C_0} \ln(F)$$
* **Monitoramento do Balanço de Massa Assimilada Acumulada:**
$$\frac{M_a}{M_m^0} = \frac{r}{1 - r}(1 - F)$$

---

## Modos de Visualização Gráfica

O simulador permite alternar em tempo real entre cinco modos de representação:

1. **Evolução Relativa Contínua ($C_l / C_0$ vs. $F$):** Curvas contínuas de enriquecimento ou empobrecimento em função da fração de líquido gerado ou remanescente.
2. **Diagrama Multielementar Relativo (**$C_l / C_0$ **vs. Elemento):** Diagramas normalizados à fonte inicial ($C_0$) para 5 estágios discretos de $F$, com cálculo automático dos fatores de fracionamento relativo $(\mathrm{La}/\mathrm{Yb})\_{l/0}$ e $(\mathrm{Eu}/\mathrm{Eu}^*)\_{l/0}$.
3. **Concentração Absoluta Contínua ($C_l$ em ppm vs. $F$):** Evolução das concentrações absolutas no líquido (ppm) a partir de um reservatório parental $C_0$.
4. **Diagrama Multielementar Normalizado ao Condrito C1 ($C_l / C_{\text{Condrito}}$):** Calcula as concentrações absolutas no líquido em ppm e normaliza ao Condrito C1 (Sun & McDonough, 1989; McDonough & Sun, 1995), eliminando o efeito Oddo-Harkins e exibindo na legenda as razões petrogenéticas reais $(\text{La/Yb})_N$ e $\text{Eu/Eu}^*$.
5. **Diagrama Multielementar em Concentração Bruta ($C_l$ em ppm vs. Elemento):** Visualização direta das abundâncias absolutas em ppm para cada estágio de $F$.

Quando apenas lantanídeos estão selecionados, o eixo horizontal ordena os elementos por número atômico ($\text{La} \to \text{Lu}$); quando LILE e HFSE estão incluídos, adota-se a sequência de incompatibilidade no manto primitivo de Sun & McDonough (1989).

---

## Banco de Dados Interno e Interoperabilidade CSV

* **Reservatórios de Referência ($C_0$ e $C_a$ em ppm):** Inclui composições pré-configuradas para **Manto Primitivo**, **Basalto N-MORB** e **Basalto OIB** (Sun & McDonough, 1989; McDonough & Sun, 1995), além de **Crosta Continental Superior** (Rudnick & Gao, 2003).
* **Base de Coeficientes de Partição ($K_d$):** Compilação didática baseada em Arth (1976), Fujimaki et al. (1984), Rollinson (1993) e Rollinson & Pease (2021), com consistência para o par geoquímico CHARAC ($\text{Y–Ho}$) e continuidade na parábola de Onuma.
* **Edição e Importação Universal de CSV:** Todas as células de $K_d$, $C_0\text{ (ppm)}$ e $C_a\text{ (ppm)}$ são editáveis na interface. O leitor de arquivos CSV detecta automaticamente tanto o padrão internacional (separador `,` e decimal `.`) quanto arquivos exportados pelo Microsoft Excel em português (separador `;` e decimal `,`).
* **Exportação Reprodutível:** Exporta matrizes de dados em `.csv` com cabeçalho completo de metadados ($X_{0,i}$, $p_i$, $D_0$, $P$, $F_{\max}$, $r$, $C_0$, $C_a$) e salva figuras prontas para publicação em `.png` (300 DPI), `.pdf` e `.svg`.

---

## Instalação e Execução

### Opção 1: Executável Standalone para Windows (Sem necessidade de instalar Python)
1. Acesse a seção [**Releases**](https://github.com/ocaminhodaspedras/Magma_Trace/releases) deste repositório.
2. Baixe o arquivo `Magma_Trace.exe` da versão mais recente e execute-o diretamente.
3. Alternativamente, solicite uma cópia para felipe.leitzke@ufrgs.br

### Opção 2: Execução via Código-Fonte (Python 3.9+)
Clone o repositório, instale as dependências listadas em `requirements.txt` e execute o script principal:

```bash
git clone [https://github.com/ocaminhodaspedras/Magma_Trace.git](https://github.com/ocaminhodaspedras/Magma_Trace.git)
cd Magma_Trace
pip install -r requirements.txt
python magma_trace.py
