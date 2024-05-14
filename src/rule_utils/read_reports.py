import json

def read_hostile_report(filepath) -> dict:
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data[0]


def get_reads_into_hostile(hostile_dict) -> int:
    return hostile_dict["reads_in"]


def get_nonhost_reads(hostile_dict) -> int:
    return hostile_dict["reads_out"]


def get_percent_host(hostile_dict) -> float:
    return hostile_dict["reads_removed_proportion"]


def main():
    print("Not yet implemented")


if __name__ == "__main__":
    main()