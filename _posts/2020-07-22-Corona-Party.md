---
layout: post
title: "Estimating the risks of partying during a pandemic"
date: 2020-07-22 10:30:00 +0100
categories: R
status: publish
published: true
# status: development
# published: false
---
 

 
*This blog post was originally published on July $22^{\text{th}}$, but was updated on August $9^{\text{th}}$ to compare the risks of partying in Amsterdam, Barcelona, and London using the most recent coronavirus case numbers.*
 
There is no doubt that, every now and then, one ought to celebrate life. This usually involves people coming together, talking, laughing, dancing, singing, shouting; simply put, it means throwing a party. With temperatures rising, summer offers all the more incentive to organize such a joyous event. Blinded by the light, it is easy to forget that we are, unfortunately, still in a pandemic. But should that really deter us?
 
Walking around central Amsterdam after sunset, it is easy to notice that not everybody holds back. Even if my Dutch was better, it would likely still be difficult to convince groups of twenty-somethings of their potential folly. Surely, they say, it is exceedingly unlikely that this little party of ours results in any virus transmission?
 
Government retorts by shifting perspective: while the chances of virus spreading at any one party may indeed be small, this does not licence throwing it. Otherwise many parties would mushroom, considerably increasing the chances of virus spread. Indeed, government stresses, this is why such parties remain *illegal*.
 
But while *if-everybody-did-what-you-did* type of arguments score high with parents, they usually do no score high with their children. So instead, in this post, we ask the question from an individual's perspective: what are the chances of getting the virus after attending this or that party? And what factors make this more or less likely?
 
As a disclaimer, I should say that I am not an epidemiologist --- who, by the way, are a [more cautious bunch](https://www.nytimes.com/interactive/2020/06/08/upshot/when-epidemiologists-will-do-everyday-things-coronavirus.html) than I or the majority of my age group --- and so my assessment of the evidence may not agree with expert opinion. With that out of the way, and without further ado, let's dive in.
 
 
# Risky business?
To get us started, let's define the *risk of a party* as the probability that somebody who is infected with the novel coronavirus and can spread it attends the gathering. The two major factors influencing this probability are the size of the party, that is, the number of people attending the gathering; and the prevalence of infectious people in the relevant population. As we will see, the latter quantity is difficult to estimate. The probability of actually getting infected by a person who has the coronavirus depends further on a number of factors; we will discuss those in a later section.
 
Let's compare the risk of partying across three wonderful European cities: Amsterdam, Barcelona, and London. From July $22^{\text{nd}}$ to August $4^{\text{th}}$, a total of $563$, $3301$, and $1101$ new infections were reported (see [here](https://www.rivm.nl/en/novel-coronavirus-covid-19/current-information), [here](https://dadescovid.cat/diari?drop_es_residencia=2&tipus=regio&id_html=ambit_2&codi=13), [here](https://coronavirus.data.gov.uk/), and the *Post Scriptum*). This results in a relative case count of $64.50$, $203.72$, and $12.54$ per $100,000$ inhabitants, respectively. While these are the numbers of *reported new infected* cases, they are not the numbers of *currently infectious* cases. How do we arrive at those?
 
 
# Estimating the true number of infectious cases
Upon infection, it usually takes a while until one can infect others, with [estimates ranging](https://theconversation.com/how-long-are-you-infectious-when-you-have-coronavirus-135295) from $1$ - $3$ days before showing symptoms. The *incubation period* is the time it takes from getting infected to showing symptoms. It lasts about $5$ days on average, with the vast majority of people showing symptoms within $12$ days (Lauer et al., [2020](https://www.acpjournals.org/doi/10.7326/M20-0504)). Yet about a third to a half of people can be infectious without showing any symptoms (Pollán et al. [2020](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)31483-5); He et al. [2020](https://www.nature.com/articles/s41591-020-0869-5)). [Estimates suggest](https://theconversation.com/how-long-are-you-infectious-when-you-have-coronavirus-135295) that one is infectious for about $8$ - $10$ days, but it can be longer.
 
These are complications, but we need to keep it simple. Currently, visitors from outside Europe must show a negative COVID-19 test or need to self-isolate for $14$ days upon arrival in most European countries (see [Austria](https://www.austria.info/en/service-and-facts/coronavirus-information), for an example). Let's take these $14$ days for simplicity, and assume conservatively that this is the time one is infectious upon getting infected. Thus, we simply take the reported number of *new infected* cases in the last two weeks as the reported number of *currently infectious* cases.[^1]
 
We have dealt with the first complication, but a second one immediately follows: how do we get from the *reported* number of infections to the *true* number of infections? One can estimate the true number of infections using models, or by empirically estimating the seroprevalence in the population, that is, the proportion of people who have developed antibodies.
 
Using the first approach, Flaxman et al. ([2020](https://www.nature.com/articles/s41586-020-2405-7)) estimate the total percentage of the population that has been infected --- the *attack rate* --- across $11$ European countries as of May $4^{\text{th}}$. The Netherlands was, unfortunately, not included in these estimates, and so we focus on Spain and the UK. For these countries the estimated true attack rates were $5.50\%$ and $5.10\%$, respectively. Given the population of these countries and the cumulative number of reported infections, we can compute the *reported* attack rate. Relating this to the estimate of the *true* attack rate gives us an indication of the extent that the reports undercount the actual infections; the code below calculates this for the three countries.
 

{% highlight r %}
library('dplyr')
library('COVID19')
 
get_undercount <- function(country, attack_rate) {
  
  # Flaxman et al. (2020) estimate the attack rate as of 4th of May
  dat <- covid19(country = country, end = '2020-05-04', verbose = FALSE)
  
  dat %>% 
    group_by(id) %>% 
    summarize(
      population = head(population, 1),
      total_cases = tail(confirmed, 1)
    ) %>% 
    mutate(
      attack_rate = attack_rate,
      reported_attack_rate = 100 * total_cases / population,
      undercount_factor = attack_rate / reported_attack_rate
    )
}
 
get_undercount(c('Spain', 'United Kingdom'), c(5.5, 5.10))
{% endhighlight %}



{% highlight text %}
## # A tibble: 2 x 6
##   id    population total_cases attack_rate reported_attack_rate undercount_factor
##   <chr>      <int>       <int>       <dbl>                <dbl>             <dbl>
## 1 ESP     46796540      218011         5.5                0.466              11.8
## 2 GBR     66460344      191843         5.1                0.289              17.7
{% endhighlight %}
 
The table above shows that cases were undercounted by a factor of about $12$ in Spain and $18$ in the UK. The Netherlands undercounted cases by a factor of about $10$ in April (Luc Coffeng, personal communication). The attack rate estimate for Spain is confirmed by a recent seroprevalence study, which finds a similarly low overall proportion of people who have developed antibodies (around $5\%$, with substantial geographical variability) in the period between April $27^{\text{th}}$ and May $11^{\text{th}}$ (Pollán et al. [2020](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)31483-5)). In another seroprevalence study, Havers et al. [2020](https://jamanetwork.com/journals/jamainternalmedicine/fullarticle/2768834) find that between March $23^{\text{rd}}$ and May $12^{\text{nd}}$, reported cases from several areas in the United States undercounted true infections by a factor between $6$ and $24$.
 
Currently, the pandemic is not as severe in Europe as it was back when the above studies were conducted. Most importantly, the testing capacity has been ramped up in most countries. For example, while the proportion of positive tests in the Netherlands and the UK were $9.40\%$ and $6.60\%$ on the $4^{\text{th}}$ of May (the end date used in the Flaxman et al. 2020 study), they currently are $1.70\%$ and $0.60\%$, using most recent data from [here](https://ourworldindata.org/coronavirus-testing) at the time of writing. Spain's coronavirus cases peaked roughly two weeks earlier than those of the Netherlands and the UK, and so by May $4^{\text{th}}$ they had a positivity rate of $2.60\%$. By July $30^{\text{th}}$ --- the date the most recent data is available at the time of writing --- their positivity rate has nearly doubled, to $5\%$.
 
Thus, while the Netherlands and the UK seem to be tracking the epidemic much more closely, which in turn implies that the factor by which they undercount the true cases is likely lower than it was previously, Spain seems to be actually doing *worse*.
 
Cases [are rising again](https://ourworldindata.org/coronavirus/country/united-kingdom?country=NLD~ESP~GBR). Let's assume therefore that the true number of *infectious* cases is $5$ times higher then the number of reported *infected* cases. For simplicity, we assume the same factor for all countries, although Spain is likely undercounting the number of true cases by a larger factor than both the Netherlands and the UK. We assume the estimated relative *true* number of *currently infectious* cases to therefore be $5 \times 64.50 = 322.50$, $5 \times 203.72 = 1018.60$, and $5 \times 12.54 = 62.70$ per $100,000$ residents in Amsterdam, Barcelona, and London, respectively. This includes asymptomatic carriers or those that are pre-symptomatic but still can spread the virus. We will assess how robust our results are against this particular correction factor later; in the next section, we estimate the risk of a party.
 
 
# Estimating the risk of a party
What are the chances that a person who attends your party has the coronavirus and is infectious? To calculate this, we assume that party guests form an independent random sample from the population. We will discuss the implications of this crude assumption later; but for now, it allows us to estimate the desired probability in a straightforward manner.
 
Take Amsterdam as an example. There were $64.50$ reported new cases per $100,000$ inhabitants between July $22^{\text{nd}}$ and August $4^{\text{th}}$. As discussed above, we take $5 \times 64.50 = 322.50$ to be the number of *true infectious cases* per $100,000$ inhabitants. Assuming that the probability of infection is the same for all citizens (more on this later), this results in $322.50 / 100,000 = 0.003225$, which gives a $0.3225\%$ or $1$ in $310$ chance that a *single* party guest has the virus and can spread it.
 
A party with just one guest would be --- *intimate*. So let's invite a few others. What are the chances that *at least one* of them can spread the virus? We compute this by first computing the complement, that is, the probability that *no* party guest is infectious.
 
The chance that any one person from Amsterdam is not infectious is $1 - 0.003225 = 0.996775$, or $99.68\%$. With our assumption of guests forming an independent random sample from the population, the probability that none of the $n$ guests can spread the virus is $0.996775^n$.
 
In our simple calculations, the chances of at least one infectious guest showing up depends only on the size of the party and the number of true infectious cases. The figure below visualizes how these two factors interact to give the risk of a party (see Lachmann & Fox, [2020](https://www.santafe.edu/research/projects/transmission-sfi-insights-covid-19), for a similar analysis regarding school reopenings).
 

 
<img src="/assets/img/2020-07-22-Corona-Party.Rmd/risk plot sensitivity-1.png" title="plot of chunk risk plot sensitivity" alt="plot of chunk risk plot sensitivity" style="display: block; margin: auto;" />
 
Let's take a moment to unpack this figure. Each coloured line represents a combination of estimated true number of infectious cases and party size that yields the same party risk. For example, attending a party of size $20$ when the true number of infectious cases per $100,000$ inhabitants is $50$ yields a party risk of about $1\%$, but so would, roughly, attending a party of size $10$ when the true relative number of infectious cases is $100$. Thus, there is a trade-off between the size of the party and the true number of infectious cases.
 
You can get a quick overview of the risks of parties of different sizes for a fixed number of true infectious cases by checking when the gray solid lines *verticallly* cross the coloured lines. Similarly, you can get a rough understanding for the risks of a party of fixed size for different numbers of true infectious cases by checking when the gray and coloured lines cross *horizontally*. The dotted vertical lines in the figure gives our previous estimate of the true number of infectious cases for London, Barcelona, and Amsterdam.
 
What's the risk of partying in those three cities? For gatherings of size $10$, $25$, and $50$, the probability that at least one guest arrives infectious is $0.63\%$, $1.56\%$, and $3.09\%$ for London. For Amsterdam, the risks are substantially higher, with $3.18\%$, $7.76\%$, and $14.91\%$. Barcelona performs worst, with staggering risks of $9.73\%$, $22.58\%$, and $40.07\%$. These numbers are sobering, and I want you to take a moment to let them sink in. We will discuss the assumptions we had to make in order to arrive at them in the next section.[^2]
 

 
Do you happen to live neither in London, Amsterdam, nor Barcelona? Regardless of your area of residence, the figure allows you to estimate the party risk; just look up the local number of new cases in the last two weeks, multiply with a correction factor (we used $5$), and --- making the assumptions we have made so far --- the plot above gives you the probability that at least one party guest will turn up infectious with the coronavirus. The assumptions we have made are very simplistic, and indeed, if you have a more elaborate way of estimating the number of currently infectious cases, then you can use that number combined with the figure to estimate the party risk.
 
Take Rome and Berlin, for example. From July $22^{\text{nd}}$ to August $4^{\text{th}}$, they had $4.82$ and $15.63$ cases per $100,000$ inhabitants, respectively (see the *Post Scriptum*). Making the same assumptions as with the other cities, and using a correction factor of $5$, the probability of having at least one infectious guest attending a party of size $25$ are $0.60\%$ for Rome and $1.93\%$ for Berlin, respectively. The table below gives the risk for parties of size $10$, $25$, $50$, and $100$ in the five European cities.
 

{% highlight text %}
##     Rome London Berlin Amsterdam Barcelona
## 10  0.24   0.63   0.78      3.18      9.73
## 25  0.60   1.56   1.93      7.76     22.58
## 50  1.20   3.09   3.83     14.91     40.07
## 100 2.38   6.08   7.52     27.60     64.08
{% endhighlight %}
 
While we have computed the party risk for a single party, this risk naturally increases when you attend multiple ones. Suppose you have been invited to parties of size $20$, $35$, and $50$ which will take place in the next month. Let's for simplicity assume that all guests are different each time. Let's further assume that the number of infectious cases stays constant over the next month. Together, these assumptions allow us to calculate the *total* party risk as the party risk of attending a single party of size $20 + 35 + 50 = 105$, which gives a considerable risk of $6.37\%$ for London, a whopping $28.76\%$ for Amsterdam, and a crippling $65.87\%$ for Barcelona. It seems that, in this case, fortune does not favour the bold.
 

 
 
## Assumptions
The analysis above is a very rough *back-of-the-envelope* calculation. We have made a number of crucial assumptions to arrive at some numbers. That's useful as a first approximation; now we have at least some intuition for the problem, and we can critically discuss the assumptions we made. Most importantly, do these assumptions lead to *overestimates* or *underestimates* of the party risk?
 
 
### Independence
First, and most critically, we have assumed that party guests are *randomly and independently* drawn from the population. It is this assumption that allowed us to compute the joint probability that none of the party guests have the virus by multiplying the probabilities of any individual being virus-free. If you have ever been to a party, you know that this is not true: instead, a considerable number of party guests usually know each other, and it is safe to say that they are similar on a range of socio-demographic variables such as age and occupation.[^3]
 
This means we are sampling not from the whole population, as our simple calculation assumes, but from some particular subpopulation that is well connected. Since the party guests likely share social circles or even households, the *effective* party size --- in terms of being relevant for virus transmission --- is smaller than the *actual* party size; this is because these individuals share the same risks. A party with $20$ married couples seems safer than a party with $40$ singles. This would suggest that we overestimate the risk of a party.
 
 
### Uniform infection probability
At the same time, however, our calculations assume that the risk of getting the coronavirus is evenly spread across the population. We used this fact when estimating the probability that any one person has the coronavirus as the total number of cases divided by the population size.
 
The probability of infection is not, howevever, evenly distributed. For example, Pollán et al. ([2020](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)31483-5)) report a seroprevalence of people aged $65$ or more of about $6\%$, while people aged between $20$ and $34$ showed a seroprevalence of $4.4\%$ between April $27^{\text{th}}$ to May $11^{\text{th}}$. These days, however, there is a substantial rise in young people who become infected, in the [United States](https://www.nytimes.com/2020/06/25/us/coronavirus-cases-young-people.html) and likely also in Europe. Because young people are [less likely to develop symptoms](https://www.axios.com/coronavirus-young-people-spread-5a0cd9e0-1b25-4c42-9ef9-da9d9ebce367.html), the virus can spread largely undetected.
 
Moreover, it seems to me that people who would join a party are in general more *adventurous*. This would increase the chances of an infectious person attending a party; thus, our calculation above may in fact underestimate the party risk.
 
At the same time, one would hope that people who show symptoms avoid parties. If all guests do so, then only pre-symptomatic or asymptomatic spread can occur, which would reduce the party risk by a half up to two thirds. On the flip side, people who show symptoms might get tested for COVID-19 and, upon receiving a negative test, consider it safe to attend a party. This might be foolish, however; recent estimates suggest that tests miss about $20\%$ infections for people who show symptoms (Kucirka et al, [2020](https://www.acpjournals.org/doi/10.7326/M20-1495); see also [here](https://www.advisory.com/daily-briefing/2020/07/06/negative-covid)). For people without symptoms, the test performs [even worse](https://www.theatlantic.com/science/archive/2020/06/how-negative-covid-19-test-can-mislead/613246/).
 
For parties taking place in summer, it is not unlikely that many guests engaged in holiday travel in the days or weeks before the date of the party. Since travel increases the chances of infection, this would further increase the chances that at least one party guest has contracted the coronavirus.
 
 
### Estimating true infections
We have assumed that the *reported* number of *new infected* cases in the last two weeks equals the number of *currently infectious* cases. This is certainly an approximation. Ideally, we would have a geographically explicit model which, at any point in time and space, provides an estimate of the number of infectious cases. To my knowledge, we are currently lacking such a model.
 
Note that, if the people who tested positive all self-isolate or, worse, end up in hospital, this clearly curbs virus spread compared to when they would continue to roam around. The former seems more likely. Moreover, these reported cases are likely not independent either, with outbreaks usually being localized. Similar to the fact that party guests know each other, the fact that reported cases cluster would lead us to overestimate the extent of virus spread at a party.
 
At the same time, in the Netherlands, for example, only those that show symptoms can get tested. Since about a third to a half are asymptomatic or pre-symptomatic in the sense that they spread the virus considerably before symptom onset, the reported number of cases likely gives an undercount of infectious people.
 
All these complications can be summarized, roughly, in the correction factor, which gives the extent to which we believe that the reported number of infected cases undercounts the true number of currently infectious cases. We have focused on a factor of $5$, but the figure above allows you to assess the sensitivity of the results to this particular choice.
 
For a very optimistic factor of $1$ --- this means that we do not undercount the true cases --- a party of size $30$ results in a risk of $0.37\%$ for London, $1.92\%$ for Amsterdam, and $5.93\%$ for Barcelona. A factor of $5$ results in risks of $1.82\%$, $9.24\%$, and $26.45\%$, respectively. A conservative estimate, using a factor of $10$, results in risks of $3.61\%$, $17.65\%$, and $46.07\%$. You can play around with these numbers yourself. Observe how they make you feel. Personally, given what we said above about the infection probability for young and adventurous people, I am inclined to err on the side of caution.
 

 
## Estimating the probability of infection
We have the defined the party risk as the probability that at least one party guest has the coronavirus and is infectious. If this person does not spread the virus to other guests, no harm is done.
 
This is exceedingly unlikely, however. The probability of getting infected is a function of the time one is exposed to the virus, and the amount of virus one is exposed to. [Estimates suggest](https://www.erinbromage.com/post/the-risks-know-them-avoid-them) that about $1,000$ SARS-CoV-2 infectious virus particles suffice for an infection. With breathing, about $20$ viral particles diffuse into the environment per minute; this increases to $200$ for speaking; coughing or sneezing can release $200,000,000$ (!) virus particles. These do not all fall to the ground, but instead can remain suspended in the air and fill the whole room; thus, physical distancing alone might [not be enough indoors](https://www.nytimes.com/2020/07/06/health/coronavirus-airborne-aerosols.html) (the extent of airborne transmission remains debated, however; see for example Klompas, Baker, & Rhee, [2020](https://jamanetwork.com/journals/jama/fullarticle/2768396)). It seems reasonable to assume that, when party guests are crowded in a room for a number of hours, many of them stand a good chance of getting infected if at least one guest is infectious. [Masks would help](https://www.erinbromage.com/post/what-s-the-deal-with-masks), of course; but how would I sip my Negroni, wearing one?
 
It is different outdoors. A Japanese study found that virus transmission inside was about $19$ times more likely than outside (Nishiura et al. [2020](https://www.medrxiv.org/content/10.1101/2020.02.28.20029272v2)). Analyzing $318$ outbreaks in China between January $4^{\text{th}}$ and February $11^{\text{th}}$, Quian et al. ([2020](https://www.medrxiv.org/content/10.1101/2020.04.04.20053058v1)) found that only a single one occurred outdoors. This suggests that parties outdoors should be much safer than parties indoors. Yet outdoor parties feature elements unlike other outdoor events; for example, there are areas --- such as bars or [public toilets](https://www.nytimes.com/2020/06/24/style/coronavirus-public-bathrooms.html) --- which could become spots for virus transmission. They usually attract more people, too. Our simple calculations suggest, with a correction factor of $5$, that the probability that at least one person out of $150$ has the coronavirus is a staggering $38.40\%$ in Amsterdam. While, in contrast to an indoor setting, the infected person is unlikely to infect the majority of the other guests, it seems likely that at least some guests will get the virus.
 
 
# To party or not to party?
If I do not care whether I get wet or not, I will never carry an umbrella, regardless of the chances of rain. Similarly, my decision to throw (or attend) a party requires not only an estimate of how likely it is that the virus spreads at the gathering; it also requires an assessment of how much I actually care.
 
As argued above, it is almost certain that the virus spreads to other guests if one guest arrives infectious. Noting that all guests are young, one might be tempted to argue that the cost of virus spread is low. In fact, people who party might even be helping --- *heroically* --- to build [herd immunity](https://fabiandablander.com/r/Covid-Exit.html)!
 
This reasoning is foolish on two grounds. First, while the proportion of infected people who die is very small for young people --- Salje et al. ([2020](https://science.sciencemag.org/content/369/6500/208)) estimate it to be $0.0045\%$ for people in their twenties and $0.015\%$ for people in their thirties --- the picture about the non-lethal, long-term effects of the novel coronavirus is only slowly becoming clear. For some people, recovery can be [lengthy](https://www.theatlantic.com/health/archive/2020/06/covid-19-coronavirus-longterm-symptoms-months/612679/) --- much longer than the two weeks we previously believed it would take. Known as "mild" cases, they [might not be so mild after all](https://www.theguardian.com/commentisfree/2020/jul/06/coronavirus-covid-19-mild-symptoms-who). Moreover, the potential [strange neurological effects](https://www.bbc.com/future/article/20200622-the-long-term-effects-of-covid-19-infection) of a coronavirus infection are becoming increasingly apparent. All told, party animals, even those guarded by their youth, might not shake it off so easily.
 
Suppose that, even after carefully considering the potential health dangers, one is still willing to take the chances. After all, it would be a *really* good party, and we young people usually eat our veggies --- especially in Amsterdam. The trouble with infectious diseases, though, is that they travel: while you might be happy to take a chance, you and the majority of party guests will probably not self-isolate after the event, right? If infections occur at the party, the virus is thus likely to subsequently spread to other, more vulnerable parts of the population.
 
So while you might remain unharmed after attending a party, others might not. Take the story of [Bob from Chicago](https://www.erinbromage.com/post/the-risks-know-them-avoid-them), summarizing an actual infection chain reported by Ghinai et al. ([2020](https://www.cdc.gov/mmwr/volumes/69/wr/mm6915e1.htm)):
 
<blockquote><p>"Bob was infected but didn't know. Bob shared a takeout meal, served from common serving dishes, with $2$ family members. The dinner lasted $3$ hours. The next day, Bob attended a funeral, hugging family members and others in attendance to express condolences. Within $4$ days, both family members who shared the meal are sick. A third family member, who hugged Bob at the funeral became sick. But Bob wasn't done. Bob attended a birthday party with $9$ other people. They hugged and shared food at the $3$ hour party. Seven of those people became ill.</p>
 
<p>But Bob's transmission chain wasn’t done. Three of the people Bob infected at the birthday went to church, where they sang, passed the tithing dish etc. Members of that church became sick. In all, Bob was directly responsible for infecting $16$ people between the ages of $5$ and $86$. Three of those $16$ died."</p></blockquote>
 
These events took place before much of the current corona measures were put in place, but the punchline remains: parties are a matter of public, not only individual health. Don't be like Bob.
 
 
---
 
I want to thank [Denny Borsboom](https://dennyborsboom.com/) and [Luc Coffeng](https://twitter.com/luc_coffeng) for helpful discussions. I also want to thank [Andrea Bacilieri](https://www.inet.ox.ac.uk/people/andrea-bacilieri/), [Denny Borsboom](https://twitter.com/BorsboomDenny), [Tom Dablander](https://www.facebook.com/nextgendoctors), and [Charlotte Tanis](https://twitter.com/CharlotteCTanis) for helpful comments on a previous version of this blog post.
 
---
 
 
## Post Scriptum
### Data
The code below gives the data used in the main text.
 

{% highlight r %}
library('httr')
library('dplyr')
library('COVID19')
 
# See https://coronavirus.data.gov.uk/developers-guide
endpoint <- paste0(
  'https://api.coronavirus.data.gov.uk/v1/data?',
  'filters=areaType=region;areaName=London&',
  'structure={"date":"date","newCases":"newCasesBySpecimenDate"}'
)
 
response <- GET(url = endpoint, timeout(10))
 
if (response$status_code >= 400) {
    err_msg <- http_status(response)
    stop(err_msg)
}
 
# Convert response from binary to JSON:
json_text <- content(response, 'text')
data <- jsonlite::fromJSON(json_text)$data
 
# From 22nd July to 4th August
london_dat <- data %>% 
  filter(
    date >= '2020-07-22', date <= '2020-08-04'
  )
 
london_total_cases <- sum(london_dat$newCases)
london_cases <- london_total_cases / 89.82000 # per 100,000 inhabitants
 
 
# From https://www.rivm.nl/en/novel-coronavirus-covid-19/current-information
amsterdam_cases <- 64.5
amsterdam_total_cases <- 563
 
 
# https://dadescovid.cat/diari?drop_es_residencia=2&tipus=regio&id_html=ambit_2&codi=13
barcelona_total_cases <- 1502 + 1745
barcelona_total_cases <- sum(
  c(235, 263, 63, 82, 302, 327, 279, 279, 353, 71, 99, 314, 307, 327)
)
barcelona_cases <- barcelona_total_cases / 16.20343
 
 
# COVID19 has data on Italian provinces
# https://www.nytimes.com/interactive/2020/world/europe/italy-coronavirus-cases.html
italy_dat <- covid19('Italy', level = 3) %>% 
  filter(
    # From 22nd July to 4th August
    date >= '2020-07-21', date <= '2020-08-04'
  )
 
rome_dat <- italy_dat %>%  filter(administrative_area_level_3 == 'Roma')
rome_total_cases <- sum(diff(rome_dat$confirmed))
rome_cases <- (rome_total_cases / rome_dat$population[1]) * 100000
 
 
# COVID19 has data on German states
# https://www.nytimes.com/interactive/2020/world/europe/germany-coronavirus-cases.html
germany_dat <- covid19('Germany', level = 2) %>% 
  filter(
    # From 22nd July to 4th August
    date >= '2020-07-21', date <= '2020-08-04'
  )
 
berlin_dat <- germany_dat %>%  filter(administrative_area_level_2 == 'Berlin')
berlin_total_cases <- sum(diff(berlin_dat$confirmed))
berlin_cases <- berlin_total_cases / 37.69495
 
 
c(rome_cases, london_cases, berlin_cases, amsterdam_cases, barcelona_cases)
{% endhighlight %}



{% highlight text %}
## [1]   4.821241  12.536183  15.625435  64.500000 203.722298
{% endhighlight %}
 
### Figure
The code below reproduces the figure in the main text.
 

{% highlight r %}
library('RColorBrewer')
 
# Probability that no guest is infectious
prob_virus_free <- function(n, true_relative_cases = 64.50 * 5) {
  prob_virus <- true_relative_cases / 100000
  (1 - prob_virus)^n
}
 
# Probability that at least one guest is infectious
party_risk <- function(n, true_relative_cases = 64.50 * 5) {
  1 - prob_virus_free(n, true_relative_cases)
}
 
# Calculates the party size that results in 'prob_virus_free'
# for a given 'true_relative_cases'
get_party_size <- function(prob_virus_free, true_relative_cases) {
  log(prob_virus_free) / log(1 - true_relative_cases / 100000)
}
 
 
plot_total_risk <- function(party_sizes, true_relative_cases, ...) {
  
  plot(
    true_relative_cases, ns, type = 'n', xaxs = 'i', yaxs = 'i', axes = FALSE,
    xlab = 'Estimated True Number of Infectious Cases per 100,000 Inhabitants',
    ylab = 'Party Size', ...
  )
  
  ticks_x <- seq(0, 1200, 100)
  minor_ticks_x <- seq(0, 1200, 50)
  
  ticks_y <- seq(0, 300, 50)
  minor_ticks_y <- seq(0, 300, 25)
  
  axis(1, at = ticks_x, cex.axis = 1.1)
  axis(2, at = ticks_y, las = 2, cex.axis = 1.1)
  rug(x = minor_ticks_x, ticksize = -0.01, side = 1)
  rug(x = minor_ticks_y, ticksize = -0.01, side = 2)
  
  abline(h = minor_ticks_y, col = 'gray86')
  abline(v = minor_ticks_x, col = 'gray86')
  
  probs_virus <- seq(0.01, 0.99, 0.01)
  party_sizes <- sapply(1 - probs_virus, get_party_size, true_relative_cases)
  
  cols <- rev(heat.colors(50))
  cols <- colorRampPalette(cols, bias = 3)(99)[-1]
  
  ix <- c(1, seq(5, 95, 5))
  show_text <- ix
  diagonal <- data.frame(
    'x' = seq(0, 1200, length.out = 300),
    'y' = seq(0, 1200, length.out = 300) * 3/12
  )
  
  for (i in ix) {
    y <- party_sizes[, i]
    line <- data.frame('x' = true_relative_cases[-1], 'y' = y[-1])
    lines(line, col = cols[i], lwd = 2)
    
    j <- reconPlots::curve_intersect(line, diagonal)
    
    if (i %in% show_text) {
      text(
        j$x + 1, j$y + 1, paste0(probs_virus[i] * 100, '%'),
        cex = 1
      )
    }
  }
}
 
ns <- seq(0, 300)
true_relative_cases <- seq(0, 1200, length.out = length(ns))
 
plot_total_risk(
  ns, true_relative_cases,
  main = 'Probability That at Least One Guest is Infectious',
  font.main = 1, cex.main = 1.75, cex.lab = 1.50
)
 
lwd <- 1.5
lines(c(london_cases * 5, london_cases * 5), c(0, 300), lty = 2)
arrows(122, 160, 68, 130, length = 0.10, lwd = lwd)
text(135, 167, 'London', cex = 1.50)
 
lines(c(amsterdam_cases * 5, amsterdam_cases * 5), c(0, 300), lty = 2)
arrows(390, 200, 330, 170, length = 0.10, lwd = lwd)
text(412, 207, 'Amsterdam', cex = 1.50)
 
lines(c(barcelona_cases * 5, barcelona_cases * 5), c(0, 500), lty = 2)
arrows(950, 167, 1012, 130, length = 0.10, lwd = lwd)
text(946, 173, 'Barcelona', cex = 1.50)
{% endhighlight %}
 
### Table
The code below reproduces the table in the main text.
 

{% highlight r %}
n <- c(10, 25, 50, 100)
 
rome_risk <- party_risk(n, 5 * rome_cases)
london_risk <- party_risk(n, 5 * london_cases)
berlin_risk <- party_risk(n, 5 * berlin_cases)
amsterdam_risk <- party_risk(n, 5 * amsterdam_cases)
barcelona_risk <- party_risk(n, 5 * barcelona_cases)
 
tab <- cbind(
  rome_risk, london_risk, berlin_risk, amsterdam_risk, barcelona_risk
)
 
colnames(tab) <- c('Rome', 'London', 'Berlin', 'Amsterdam', 'Barcelona')
rownames(tab) <- n
 
round(tab * 100, 2)
{% endhighlight %}



{% highlight text %}
##     Rome London Berlin Amsterdam Barcelona
## 10  0.24   0.63   0.78      3.18      9.73
## 25  0.60   1.56   1.93      7.76     22.58
## 50  1.20   3.09   3.83     14.91     40.07
## 100 2.38   6.08   7.52     27.60     64.08
{% endhighlight %}
 
---
 
## References
- Flaxman, Mishra, Gandy et al. ([2020](https://www.nature.com/articles/s41586-020-2405-7)). Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe. *Nature*, 3164.
- Ghinai, I., Woods, S., Ritger, K. A., McPherson, T. D., Black, S. R., Sparrow, L., ... & Arwady, M. A. ([2020](https://www.cdc.gov/mmwr/volumes/69/wr/mm6915e1.htm)). Community Transmission of SARS-CoV-2 at Two Family Gatherings-Chicago, Illinois, February-March 2020. *MMWR. Morbidity and mortality weekly report, 69*(15), 446.
- Havers, F. P., Reed, C., Lim, T. W., Montgomery, J. M., Klena, J. D., Hall, A. J., ... & Krapiunaya, I. ([2020](https://jamanetwork.com/journals/jamainternalmedicine/fullarticle/2768834)). Seroprevalence of Antibodies to SARS-CoV-2 in Six Sites in the United States, March 23-May 3, 2020. *JAMA Internal Medicine*.
- He, X., Lau, E. H., Wu, P., Deng, X., Wang, J., Hao, X., ... & Mo, X. ([2020](https://www.nature.com/articles/s41591-020-0869-5)). Temporal dynamics in viral shedding and transmissibility of COVID-19. *Nature Medicine, 26*(5), 672-675.
- Klompas, M., Baker, M. A., & Rhee, C. ([2020](https://jamanetwork.com/journals/jama/fullarticle/2768396)). Airborne Transmission of SARS-CoV-2: Theoretical Considerations and Available Evidence. *JAMA*.
- Kucirka, L. M., Lauer, S. A., Laeyendecker, O., Boon, D., & Lessler, J. ([2020](https://www.acpjournals.org/doi/full/10.7326/M20-1495)). Variation in false-negative rate of reverse transcriptase polymerase chain reaction–based SARS-CoV-2 tests by time since exposure. *Annals of Internal Medicine*.
- Lachmann, M., & Fox, S. ([2020](https://sfi-edu.s3.amazonaws.com/sfi-edu/production/uploads/ckeditor/2020/07/07/t-034-lachmann.pdf)). When thinking about reopening schools, an important factor to consider is the rate of community transmission. *Santa Fe Institute Transmission*.
- Lauer, S. A., Grantz, K. H., Bi, Q., Jones, F. K., Zheng, Q., Meredith, H. R., ... & Lessler, J. ([2020](https://www.acpjournals.org/doi/10.7326/M20-0504)). The incubation period of coronavirus disease 2019 (COVID-19) from publicly reported confirmed cases: estimation and application. *Annals of Internal Medicine, 172*(9), 577-582.
- Morawska, L., & Cao, J. ([2020](https://www.sciencedirect.com/science/article/pii/S016041202031254X)). Airborne transmission of SARS-CoV-2: The world should face the reality. *Environment International*, 105730.
- Nishiura, H., Oshitani, H., Kobayashi, T., Saito, T., Sunagawa, T., Matsui, T., ... & Suzuki, M. ([2020](https://www.medrxiv.org/content/10.1101/2020.02.28.20029272v2)). Closed environments facilitate secondary transmission of coronavirus disease 2019 (COVID-19). *medRxiv*.
- Pollán, M., Pérez-Gómez, B., Pastor-Barriuso, R., Oteo, J., Hernán, M. A., Pérez-Olmeda, M., ... & Molina, M. ([2020](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)31483-5/fulltext)). Prevalence of SARS-CoV-2 in Spain (ENE-COVID): a nationwide, population-based seroepidemiological study. *The Lancet*.
- Salje, H., Kiem, C. T., Lefrancq, N., Courtejoie, N., Bosetti, P., Paireau, J., ... & Le Strat, Y. ([2020](https://science.sciencemag.org/content/369/6500/208)). Estimating the burden of SARS-CoV-2 in France. *Science*.
- Qian, H., Miao, T., Li, L. I. U., Zheng, X., Luo, D., & Li, Y. ([2020](https://www.medrxiv.org/content/10.1101/2020.04.04.20053058v1)). Indoor transmission of SARS-CoV-2. *medRxiv*.
 
---
 
## Footnotes
[^1]: Reported deaths are more reliable than reported cases because deaths must always be reported. This is why, for example, Flaxman et al. ([2020](https://www.nature.com/articles/s41586-020-2405-7)) use deaths to estimate the actual proportion of infections. There are issues with reported deaths, too, however, and I discuss some of them [here](https://scienceversuscorona.com/visualising-the-covid-19-pandemic/).
[^2]: Let me note that RIVM --- the Dutch National Institute for Public Health and the Environment --- has their own estimate of the number of currently infectious cases. On August $4^{\text{th}}$, their [dashboard](https://coronadashboard.rijksoverheid.nl/) showed an estimate of $94.30$ infectious cases per $100,000$ inhabitants. This number is larger than $60.40$, the number of reported number cases per $100,000$ in Amsterdam between July $22^{\text{th}}$ and August $4^{\text{th}}$. In the terms of our calculations, their model applies a correction factor of $94.30 / 60.40 = 1.56$. RIVM is therefore slightly more optimistic than I am; for parties of size $10$, $25$, and $50$, their estimates of the probability that at least one guest is infectious --- assuming guests form a random sample from the population --- are $0.94\%$, $2.33\%$, and $4.61\%$, respectively. How does RIVM arrive at their estimate of the number of infectious cases? We currently do not know. [Their weekly report](https://www.rivm.nl/documenten/wekelijkse-update-epidemiologische-situatie-covid-19-in-nederland) (Section 9.1) devotes only two small paragraphs to it, saying that the method is "still under development".
[^3]: Once the pandemic is over, inviting a random sample from the population should definitely become a thing. Bursting bubbles, one party at a time!
