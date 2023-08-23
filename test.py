import pandas as pd

data = pd.read_csv("/home/answain/alice/CircularCrossSections/FineTuningTest/cylinder_data.csv")
print(data.head())

categories = ['positiveZ','negativeZ','mainBarrel']


positiveZ_data = data[data['DetectorPart'] == categories[0]]
negativeZ_data = data[data['DetectorPart'] == categories[1]]
main_barrel_data = data[data['DetectorPart'] == categories[2]]

positiveZ_data_sorted = positiveZ_data.sort_values(by='Zmax',ascending=False)
negativeZ_data_sorted = negativeZ_data.sort_values(by='Zmin',ascending=True)


print(positiveZ_data_sorted.head())
print(negativeZ_data_sorted.head())

