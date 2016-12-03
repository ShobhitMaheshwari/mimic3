
# coding: utf-8

# In[28]:

class Feature(object):
    def __init__(self, path="~/mimic3/mimic3/demo/"):
        self.path = path
    
    def get_items(self):
        HR = [220045 ,211] #[Heart Rate, Heart Rate]
        NIBP_S = [220179] #[Non Invasive Blood Pressure systolic]
        NIBP_D = [220180] #[Non Invasive Blood Pressure diastolic]
        CAP_Refill = [115, 3348, 8377, 223951, 224308] #[Capillary Refill [Right], Capillary Refill, Capillary Refill [Left], Capillary Refill R, Capillary Refill L]
        Glucose = [227015] #[Glucose_ApacheIV]
        pH = [227037] #[PH_ApacheIV]
        Temp = [227054] #[TemperatureF_ApacheIV]
        UrineScore = [227059] #[UrineScore_ApacheIV]
        Oxygen = [227035] #[OxygenScore_ApacheIV]
        RespiratoryRate = [224688] #[Respiratory Rate (Set)]
        GCS = [198] #[GCS Total]
        FiO2 = [1040] #[BIpap FIO2]
        ETCO2 = [1817] #[ETCO2]
#         items = [HR, NIBP_S, NIBP_D, CAP_Refill, Glucose, pH, Temp, UrineScore, Oxygen, RespiratoryRate, GCS, FiO2, ETCO2]
        return {
            'HR': HR,
            'NIBP_S': NIBP_S,
            'NIBP_D': NIBP_D,
            'CAP_Refill': CAP_Refill,
            'Glucose': Glucose,
            'pH': pH,
            'Temp': Temp,
            'UrineScore': UrineScore,
            'Oxygen': Oxygen,
            'RespiratoryRate': RespiratoryRate,
            'GCS': GCS,
            'FiO2': FiO2,
            'ETCO2': ETCO2
        }
    
    def get_relative_time(self, admit_time, chart_time):
        from datetime import datetime
        # Accepts string in "2144-07-24 09:00:00" format and returns difference in seconds
        adm_time = datetime.strptime(admit_time, '%Y-%m-%d %H:%M:%S')
        char_time = datetime.strptime(chart_time, '%Y-%m-%d %H:%M:%S')
        del_time = char_time - adm_time
        rel_time = del_time.days*86400 + del_time.seconds
        return rel_time
    
    def get_patient_last_recorded_time(self, df_patient, admit_time):
        max_time = 0
        for index, row in df_patient.iterrows():
            temp = self.get_relative_time(admit_time, row['CHARTTIME'])
            if (temp > max_time):
                max_time = temp
        return max_time
    
    def get_relative_time_per_admid(self, df, admit_time):
        rlt = []
        max_time = 0
        for cTime in df['CHARTTIME']:
            temp = self.get_relative_time(admit_time, cTime)
            rlt.append(temp)
            if (temp > max_time):
                max_time = temp
            
        df['RELTIME'] = rlt
        return (df, max_time)
    
    def get_patient_event(self, df, item_ids, from_time_seconds, to_time_seconds, old_val):
        import math
        
        m = df[df.ITEMID.isin(item_ids)]
        l = m[(m['RELTIME'] >= from_time_seconds) & (m['RELTIME'] < to_time_seconds)]
        j = l['VALUE'].astype('float64').mean()
        if math.isnan(j): 
            return [old_val, 0]
        else:
            return [j, 1]
    
    def get_features(self):
        import pandas as pd
        admission_df = pd.read_csv(self.path+'ADMISSIONS.csv', sep=',',header=0)
        chartevents_df = pd.read_csv(self.path+'CHARTEVENTS.csv')
        #hardcoded
        admission_df = admission_df[admission_df['SUBJECT_ID'] == 10006]
        chartevents_df = chartevents_df[chartevents_df['SUBJECT_ID'] == 10006]
        
        patients = admission_df.groupby('SUBJECT_ID')
        features = []
        for subject_id in patients.groups:
            patient = patients.get_group(subject_id)#patient is a dataframe for one patient

            for index, patient_1 in patient.iterrows():
                (chartevents_perAdm_df, max_time) = self.get_relative_time_per_admid(
                    chartevents_df[chartevents_df['HADM_ID'] == patient_1['HADM_ID']],
                    patient_1['ADMITTIME'])# max_time is in seconds
                
                feature_patient = []
                prev_val = {}
                for it in range(0, max_time, 60*10):
                    feature_patient_time = []
                    items = self.get_items()
                    for item_name in items:
                        prev_val[item_name] = self.get_patient_event(chartevents_perAdm_df, 
                                                items[item_name], it, it+(60*10), 
                                                           prev_val[item_name] if item_name in prev_val else 0)
                        feature_patient_time.extend(prev_val[item_name])
                    feature_patient.append(feature_patient_time)
                features.append(feature_patient)
        return features


# In[29]:

feature = Feature()
print(feature.get_features())


# In[ ]:



