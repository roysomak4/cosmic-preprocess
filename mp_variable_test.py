import multiprocessing


def main():
    # manager = multiprocessing.Manager()
    # final_list = manager.list()
    final_list = []
    
    input_list_one = ['one', 'two', 'three', 'four', 'five']
    input_list_two = ['six', 'seven', 'eight', 'nine', 'ten']


    process1 = multiprocessing.Process(target=worker, args=[input_list_one, final_list])
    process2 = multiprocessing.Process(target=worker, args=[input_list_two, final_list])

    process1.start()
    process2.start()
    process1.join()
    process2.join()

    print(final_list)


def worker(data, target_list):
    for item in data:
        target_list.append(item)


if __name__ == '__main__':
    main()
