
def merge(intervals) -> list:
    """
    merge intervals

    Parameters:
    intervals: list includes ranges like [start, end]

    Returns:
    list: merged intervals 

    Example:
    >>> merge([[1,3],[2,6],[8,10],[15,18]])
    [[1,6],[8,10],[15,18]]
    """

    intervals_sorted = sorted(intervals, key=lambda x : x[0])
    # print(intervals_sorted)
    result = []
    for interval in intervals_sorted:
       
        if result and result[-1][1] >= interval[0]:
           
            result[-1][1] = max(result[-1][1], interval[1])
            # print(result)
        else:
            result.append(interval)
    # print(result)
    return result



def intervalIntersection(A, B) -> list:

    """
    calculate intersection of intervals

    Parameters:
    A and B: two intervals list

    Returns:
    list: intersection intervals 

    Example:
    >>> intervalIntersection(A = [[0,2],[5,10],[13,23],[24,25]], B = [[1,5],[8,12],[15,24],[25,26]])
    [[1,2],[5,5],[8,10],[15,23],[24,24],[25,25]]
    """

    i, j = 0, 0 # 双指针
    res = []
    A = sorted(A, key=lambda x : x[0])
    B = sorted(B, key=lambda x : x[0])
    while i < len(A) and j < len(B):
        a1, a2 = A[i][0], A[i][1]
        b1, b2 = B[j][0], B[j][1]
        # 两个区间存在交集
        if b2 >= a1 and a2 >= b1:
            # 计算出交集，加入 res
            res.append([max(a1, b1), min(a2, b2)])
        # 指针前进
        if b2 < a2: j += 1
        else:       i += 1

    # for i in res:
    #     if (i[1]- i[0]) <= 0:
    #         print(i)
    
    return res





def calculate_genebody(file_readlines: list) -> int:

    gene_interval_dict = {}
    
    for i in file_readlines:
        if i[2] not in gene_interval_dict:
            gene_interval_dict[i[2]] = [[int(i[4]), int(i[5])]] 
        else:
            gene_interval_dict[i[2]].append([int(i[4]), int(i[5])])


    sum = 0
    for _ in gene_interval_dict.values():
        
    
    # print(merge(intervals))

        a = merge(intervals=_)
        for i in a:
            sum += (i[1] - i[0] + 1)

    # print(sum)
    return sum


def parse_exons_to_list(i: list) -> list:
    _ = []

    exon_start_string = i[-2].split(",")[0:-1]
    exon_end_string = i[-1].split(",")[0:-1]
    for j in range(len(exon_end_string)):
        interval_of_the_exon = [int(exon_start_string[j]), int(exon_end_string[j])]
        _.append(interval_of_the_exon)
        
    return _



def calculate_exons(file_readlines: list) -> int:

    exon_interval_dict = {}
    
    for i in file_readlines:
        if i[2] not in exon_interval_dict:
            
            _ = []
            exon_start_string = i[-2].split(",")[0:-1]
            exon_end_string = i[-1].split(",")[0:-1]
            for j in range(len(exon_end_string)):
                interval_of_the_exon = [int(exon_start_string[j]), int(exon_end_string[j])]
                _.append(interval_of_the_exon)
                exon_interval_dict[i[2]] = _
        else:

            exon_start_string = i[-2].split(",")[0:-1]
            exon_end_string = i[-1].split(",")[0:-1]
            for j in range(len(exon_end_string)):
                interval_of_the_exon = [int(exon_start_string[j]), int(exon_end_string[j])]
            # _.append(interval_of_the_exon)
                exon_interval_dict[i[2]].append(interval_of_the_exon)

    sum = 0
    for _ in exon_interval_dict.values():
        
    
    # print(merge(intervals))

        a = merge(intervals=_)
        for i in a:
            sum += (i[1] - i[0] + 1)

    # print(sum)
    return sum



def UTR_region_5(file_readlines) -> int:

    UTR_5_dictionary = {}
    for i in file_readlines:
        if i[3] == "+" and i[6] != i[7]:
            cds_5_UTR_region_of_gene = [[int(i[4]), int(i[6])-1]]
            exon_and_CDS_5_UTR_intersection = intervalIntersection(parse_exons_to_list(i=i), cds_5_UTR_region_of_gene)


            if i[2] not in UTR_5_dictionary:
                UTR_5_dictionary[i[2]] = exon_and_CDS_5_UTR_intersection
            else:
                for k in exon_and_CDS_5_UTR_intersection:
                    UTR_5_dictionary[i[2]].append(k)
        
        if i[3] == "-" and i[6] != i[7]:
            
            cds_5_UTR_region_of_gene = [[int(i[7])+1, int(i[5])]]
            exon_and_CDS_5_UTR_intersection = intervalIntersection(parse_exons_to_list(i=i), cds_5_UTR_region_of_gene)

            if i[2] not in UTR_5_dictionary:
                UTR_5_dictionary[i[2]] = exon_and_CDS_5_UTR_intersection
            else:
                for k in exon_and_CDS_5_UTR_intersection:
                    UTR_5_dictionary[i[2]].append(k)
    
  
    sum = 0

    for _ in UTR_5_dictionary.values():
    # print(merge(intervals))
        a = merge(intervals=_)
        # a = _
        for i in a:
            sum += (i[1] - i[0] + 1)



    return sum



def non_coding(file_readlines) -> int:

    
    for i in file_readlines:
        if i[6] == i[7]:
            gene_interval_dict = {}

            if i[2] not in gene_interval_dict:
                gene_interval_dict[i[2]] = [[int(i[4]), int(i[5])]] 
            else:
                gene_interval_dict[i[2]].append([int(i[4]), int(i[5])])

    sum = 0
    for _ in gene_interval_dict.values():
    # print(merge(intervals))

        a = merge(intervals=_)
        for i in a:
            
            sum += (i[1] - i[0] + 1)

    # print(sum)
    return sum




def UTR_region_3(file_readlines) -> int:

    UTR_5_dictionary = {}
    for i in file_readlines:
        if i[3] == "+" and int(i[6]) != int(i[7]):
            cds_5_UTR_region_of_gene = [[int(i[7])+1, int(i[5])]]
            exon_and_CDS_5_UTR_intersection = intervalIntersection(cds_5_UTR_region_of_gene, parse_exons_to_list(i=i))


            if i[2] not in UTR_5_dictionary:
                UTR_5_dictionary[i[2]] = exon_and_CDS_5_UTR_intersection
            else:
                for k in exon_and_CDS_5_UTR_intersection:
                    UTR_5_dictionary[i[2]].append(k)
        
        if i[3] == "-" and int(i[6]) != int(i[7]):
            
            cds_5_UTR_region_of_gene = [[int(i[4]), int(i[6])-1]]
            exon_and_CDS_5_UTR_intersection = intervalIntersection(cds_5_UTR_region_of_gene, parse_exons_to_list(i=i))

            if i[2] not in UTR_5_dictionary:
                UTR_5_dictionary[i[2]] = exon_and_CDS_5_UTR_intersection
            else:
                for k in exon_and_CDS_5_UTR_intersection:
                    UTR_5_dictionary[i[2]].append(k)
    


    sum = 0
    for _ in UTR_5_dictionary.values():
    # print(merge(intervals))
        a = merge(intervals=_)
        # a = _
        for i in a:
            sum += (i[1] - i[0] + 1)

    # print(sum)
    return sum


if '__main__' == __name__:

    with open("C://Users/wangt/Desktop/refFlat.txt", "r") as ref:
        first_line = ref.readlines()
        first_line = [i.strip().split("\t") for i in first_line]

    
    
    print(calculate_exons(file_readlines=first_line))
    print(calculate_genebody(first_line))


    print(UTR_region_5(first_line))
    print(UTR_region_3(file_readlines=first_line))