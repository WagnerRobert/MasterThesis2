__author__ = 'delur'

def checkOrder(svm_dict, result):
#    print svm_dict
    sum_list = []

    for protein in svm_dict:
        total = 0.0
        for kmer, value in svm_dict[protein]:
            total += value
        sum_list.append( (protein, total))


    sum_list = sorted(sum_list, key=lambda  x:x[1])

    for name, total in sum_list:
        print result[name][0] + "\t" + str(total) +"\t" + name