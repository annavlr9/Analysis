```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tensorflow as tf
import keras

from numpy import loadtxt
from keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.preprocessing import StandardScaler, LabelEncoder
from pandas.plotting import scatter_matrix

from sklearn.metrics import auc, roc_curve, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix, classification_report, cohen_kappa_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, KFold, cross_val_score, StratifiedKFold
from sklearn.datasets import make_classification
```


```python
expression = pd.read_csv('C:/Users/annav/Documents/ML_project/data/final_data/new_alignment/annotate.txt', sep = '\t')  # Change to delimiter=' ' if space-delimited
expression.head()
```


```python
meth_CGI = pd.read_csv('C:/Users/annav/Documents/ML_project/data/final_data/new_alignment/CG_5hmC_ATAC_CpGI_final.txt', sep = '\t') 
meth_CGI.columns = ['chrom', 'start', 'end', 'Level_5hmC', 'annotation', 'ensembl', 'SYMBOL']
meth_CGI.head()
```


```python
# Merge df2 into df1 by the common column 'id'
data_promoters_CGI = pd.merge(meth_CGI, expression, on='ensembl', how='inner') 
data_promoters_CGI = data_promoters_CGI[['Level_5hmC', 'exp1', 'exp2']]
data_promoters_CGI.head()
```


```python
data_promoters_CGI.shape
```


```python
data_promoters_CGI['Category'] = data_promoters_CGI['Level_5hmC'].apply(lambda x: 'Low' if x < 0.065 else 'Medium' if x < 0.14 else 'High')
data_promoters_CGI.head()
```


```python
data_promoters_CGI.shape
```


```python
# Count the occurrences of each category
category_counts = data_promoters_CGI['Category'].value_counts()
print(category_counts)
```


```python
# PIE CHART

# Define the data
labels = ['Low', 'Medium', 'High']
sizes = [2334074, 40574, 17343]
colors = ['#1ABED2', '#2CA02C', '#FE6A0C']  # Define custom colors (light red, light blue, light green)

# Create the pie chart
plt.figure(figsize=(7, 7))  # Optional: Adjust the size of the plot
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, 
        colors=colors, textprops={'fontsize': 14, 'fontweight': 'bold'})  # Custom colors and text styling

# Add a title
plt.title('Promters CGI', fontsize=16)

# Show the pie chart
plt.show()
```


```python
# BARPLOT

Data = ['Low', 'Medium', 'High']
Num = [2334074, 40574, 17343]

# Define a list of colors for the bars
colors = ['#ff9999', '#66b3ff', '#c2f0c2']

# Create the bar plot with different colors
plt.figure(figsize=(10, 6))  # Adjust figure size
bars = plt.bar(Data, Num, color=colors)

# Add values on top of each bar
for i, value in enumerate(Num):
    plt.text(i, value + max(Num) * 0.01,  # Offset the text slightly above the bar
             str(value), ha='center', va='bottom', fontsize=10)

# Add labels and title
plt.title('', fontsize=14)
plt.xlabel('Categories', fontsize=12)
plt.ylabel('Number of CGs', fontsize=12)
plt.xticks(rotation=45, fontsize=10, ha='right')

# Show the plot
plt.tight_layout()  # Adjust layout to prevent clipping
plt.show()
```


```python
# Check for missing values
print('Missing values per column:\n', data_promoters_CGI.isnull().sum())
```


```python
## Convert the Category column to binary METHOD1

# Sample DataFrame
df_promoters_CGI = pd.DataFrame(data_promoters_CGI)

# Initialize LabelEncoder
label_encoder_promoters_CGI = LabelEncoder()

# Fit and transform the 'Category' column
df_promoters_CGI['Category_Encoded'] = label_encoder_promoters_CGI.fit_transform(df_promoters_CGI['Category'])
print(df_promoters_CGI)
df_promoters_CGI.head(20)
```


```python
## Separate features and target variable
array = df_promoters_CGI.values #The `.values` attribute returns the data as a NumPy array
print(array[:5])
```


```python
X = array[:,1:2]
y = array[:,3]
```


```python
print(X[:5])
```


```python
print(y[:5])
```


```python
# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 42)
```


```python
print(X_train[:5])
print(y_train[:5])

print(X_test[:5])
print(y_test[:5])
```


```python
# BOXPLOT INTEGRATING ALL MODELS

from sklearn.model_selection import StratifiedKFold
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, Sequential
import matplotlib.pyplot as plt


## Separate features and target variable
array = df_promoters_CGI.values #The `.values` attribute returns the data as a NumPy array

X = array[:,1:2]
y = array[:,3]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 42)




models = []
models.append(('LR', LogisticRegression(solver = 'lbfgs', multi_class = 'multinomial')))
models.append(('KNN', KNeighborsClassifier()))
models.append(('RF', RandomForestClassifier()))
models.append(('DT', DecisionTreeClassifier()))
models.append(('NB', GaussianNB()))
#models.append(('SVM', SVC(gamma = 'auto')))


label_encoder = LabelEncoder()
y_transform = label_encoder.fit_transform(y_train)

# Evaluate each model and compare the results
results = []
names = []

for name, model in models:
    kfold = StratifiedKFold(n_splits = 10, random_state = 1, shuffle = True)
    if name == "RF":
        cv_results = cross_val_score(model, X_train, y_transform, cv = kfold, scoring = 'accuracy')
    else:
        cv_results = cross_val_score(model, X_train, y_train, cv = kfold, scoring = 'accuracy')
                                    
    #cv_results = cross_val_score(model, X_train, y_train, cv=kfold, scoring='accuracy')
    results.append(cv_results)
    names.append(name)
    print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))
```


```python
# DL model
X_DP = array[:, 1:2]
y_DP = array[:, 4]

# Convert data to float32
X_DP = np.asarray(X_DP).astype(np.float32)
y_DP = np.asarray(y_DP).astype(np.int32)  # Ensure y is integer for sparse_categorical_crossentropy

# Initialize cross-validation
kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)
dl_accuracies = []

# Cross-validate DL model
for train_idx, val_idx in kfold.split(X_DP, y_DP):
    X_train, X_val = X_DP[train_idx], X_DP[val_idx]
    y_train, y_val = y_DP[train_idx], y_DP[val_idx]
    
    # Scale the data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    
    # Build the model
    model = Sequential([
        layers.Input(shape=(1,)),
        layers.Dense(10),
        layers.ReLU(),
        layers.Dense(10),
        layers.ReLU(),
        layers.Dense(10),
        layers.ReLU(),
        layers.Dense(3),
        layers.Softmax()
    ])
    
    model.compile(loss='sparse_categorical_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])
    
    # Train the model
    model.fit(X_train_scaled, y_train, epochs=50, verbose=0)
    
    # Evaluate the model on the validation fold
    _, accuracy = model.evaluate(X_val_scaled, y_val, verbose=0)
    dl_accuracies.append(accuracy) 

# Add DL results to the comparison
results.append(dl_accuracies)
names.append('DL')
```


```python
# Print DL results
print('DL: %f (%f)' % (np.mean(dl_accuracies), np.std(dl_accuracies)))
```


```python
# Plot the boxplot with all the models
plt.boxplot(results, labels=names)
plt.title('Algorithm Comparison (3 HD)')
plt.ylabel('Accuracy')
plt.show()
```


```python
# AUC PLOT INTEGRATING ALL MODELS

from sklearn.model_selection import StratifiedKFold
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, Sequential
import matplotlib.pyplot as plt
from sklearn.base import clone
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification


## Separate features and target variable
array = df_promoters_CGI.values #The `.values` attribute returns the data as a NumPy array

X = array[:,1:2]
y = array[:,3]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 42)

models = []
models.append(('LR', LogisticRegression(solver = 'lbfgs', multi_class = 'multinomial')))
models.append(('KNN', KNeighborsClassifier()))
models.append(('RF', RandomForestClassifier()))
models.append(('DT', DecisionTreeClassifier()))
models.append(('NB', GaussianNB()))


label_encoder = LabelEncoder()
y_transform = label_encoder.fit_transform(y_train)

# Evaluate each model and compare the results
results = []
names = []

for name, model in models:
    kfold = StratifiedKFold(n_splits = 10, random_state = 1, shuffle = True)
    if name == "RF":
        cv_results = cross_val_score(model, X_train, y_transform, cv = kfold, scoring = 'accuracy')
    else:
        cv_results = cross_val_score(model, X_train, y_train, cv = kfold, scoring = 'accuracy')
                                    
    #cv_results = cross_val_score(model, X_train, y_train, cv=kfold, scoring='accuracy')
    results.append(cv_results)
    names.append(name)
    print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))









# Prepare data
X_DP = array[:, 1:2]
y_DP = array[:, 4]


# Convert data to float32
X_DP = np.asarray(X_DP).astype(np.float32)
y_DP = np.asarray(y_DP).astype(np.int32)

# Binary or multi-class classification
n_classes = len(np.unique(y_DP))
y_binarized = label_binarize(y_DP, classes=np.arange(n_classes))

# Initialize cross-validation
kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)

# To store ROC data
roc_data = {}

# Machine Learning Models AUC and ROC
for name, model in models:
    model = clone(model)  # Clone the model to avoid side effects
    y_pred_prob = cross_val_predict(model, X_DP, y_DP, cv=kfold, method='predict_proba')
    auc = roc_auc_score(y_binarized, y_pred_prob, multi_class='ovr')
    print(f"{name}: AUC = {auc:.4f}")
    
    # Calculate ROC Curve
    fpr, tpr, _ = roc_curve(y_binarized.ravel(), y_pred_prob.ravel())
    roc_data[name] = (fpr, tpr, auc)

# Deep Learning Model AUC and ROC
dl_tprs = []
dl_aucs = []
mean_fpr = np.linspace(0, 1, 100)

for train_idx, val_idx in kfold.split(X_DP, y_DP):
    X_train, X_val = X_DP[train_idx], X_DP[val_idx]
    y_train, y_val = y_DP[train_idx], y_DP[val_idx]

    # Scale the data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)

    # Build the DL model
    model = Sequential([
        layers.Input(shape=(1,)),
        layers.Dense(10),
        layers.ReLU(),
        layers.Dense(10),
        layers.ReLU(),
        layers.Dense(n_classes),
        layers.Softmax()
    ])

    model.compile(loss='sparse_categorical_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])

    # Train the model
    model.fit(X_train_scaled, y_train, epochs=50, verbose=0)

    # Predict probabilities
    y_pred_prob = model.predict(X_val_scaled)

    # Compute AUC
    auc = roc_auc_score(label_binarize(y_val, classes=np.arange(n_classes)), y_pred_prob, multi_class='ovr')
    dl_aucs.append(auc)

    # Compute ROC curve
    fpr, tpr, _ = roc_curve(label_binarize(y_val, classes=np.arange(n_classes)).ravel(), y_pred_prob.ravel())
    dl_tprs.append(np.interp(mean_fpr, fpr, tpr))
    dl_tprs[-1][0] = 0.0

# Aggregate DL ROC Curve
mean_tpr = np.mean(dl_tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = np.mean(dl_aucs)
roc_data['DL'] = (mean_fpr, mean_tpr, mean_auc)

# Plot AUC-ROC Curves for All Models
plt.figure(figsize=(10, 8))
plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random Chance')

for name, (fpr, tpr, auc) in roc_data.items():
    plt.plot(fpr, tpr, label=f"{name} (AUC = {auc:.4f})")

plt.title('ROC Curves for All Models')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.grid()
plt.show()
```
