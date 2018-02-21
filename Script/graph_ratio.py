import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

np.random.seed(0)
# sns.set(style="white", color_codes=True)
sns.set(color_codes=True)
data = pd.read_csv("../Result/Align/test.txt",header=0)

observed_ratio = data['observed_ratio'].values
theorical_ratio = data['theorical_ratio'].values


observed_ratio=np.multiply(observed_ratio,100)
theorical_ratio=np.multiply(theorical_ratio,100)


#simple graph
#tips = sns.load_dataset("../Result/Align/text.txt")
g = sns.JointGrid(x=observed_ratio, y=theorical_ratio, data=data)
g = g.plot(sns.regplot, sns.distplot)
g = g.annotate(stats.pearsonr)
fig = g.fig
# fig = sns.get_figure()
fig.savefig("test.png")
exit()

# sns.set(style="ticks")
# dataset = data[['subsample', 'observed_ratio',  'therorical_ratio']]
# dataset.loc[:,'observed_ratio'] = dataset['observed_ratio']*100
# dataset.loc[:,'therorical_ratio'] = dataset['therorical_ratio']*100
# print(dataset)
# sns.lmplot(x="observed_ratio", y="therorical_ratio", col="subsample", hue="dataset", data=dataset.values)
#            # col_wrap=2, ci=None, palette="muted", size=3,
#            # scatter_kws={"s": 50, "alpha": 1}


# import seaborn as sns
# #sns.set(style="ticks")

# # Load the example dataset for Anscombe's quartet
# df = sns.load_dataset("anscombe")
dataset = data[['subsample', 'observed_ratio',  'therorical_ratio']]
dataset.loc[:,'observed_ratio'] = dataset['observed_ratio']*100
dataset.loc[:,'therorical_ratio'] = dataset['therorical_ratio']*100
# print(df)
# # Show the results of a linear regression within each dataset

# g = sns.lmplot(x="observed_ratio", y="therorical_ratio", hue="subsample", data=dataset,
#            markers=["o", "x", "s"], palette="Set1");


# # g = sns.lmplot(x="observed_ratio", y="therorical_ratio", col="subsample", hue="subsample", data=dataset,
# # 	col_wrap=2, ci=None, palette="muted",
# # 	scatter_kws={"s": 50, "alpha": 1}, kind="reg")
# # g = g.plot(sns.regplot, sns.distplot)
# # g = g.annotate(stats.pearsonr)
# fig = g.fig
# fig.savefig("test.png")

sns.set(style="ticks")
d = sns.boxplot(x="observed_ratio", y="therorical_ratio", hue="subsample", data=dataset, palette="PRGn")
sns.despine(offset=10, trim=True)
fig = d.fig
fig.savefig("test1.png")
