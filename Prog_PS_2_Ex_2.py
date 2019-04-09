class Account:

    def __init__(self, name, account_number, initial_amount, transactions):
        self.name = name
        self.no = account_number
        self.balance = initial_amount
        self.transactions = transactions

    def deposit(self, amount):
        self.balance += amount
        self.transactions += 1

    def withdraw(self, amount):
        self.balance -= amount
        self.transactions += 1

    def dump(self):
        s = "Account %s belongs to %s. Its current balance is %s CHF, after %s transactions this \
period" % (self.no, self.name, self.balance, self.transactions)
        print(s)


a1 = Account("Jon Olsson", "19371554951", 20000, 0)
a2 = Account("Liz Olsson", "19371564761", 20000, 0)
a3 = Account("Jeremie Heitz", "19371574871", 60000, 0)
a1.deposit(1000)
a1.withdraw(4000)
a2.withdraw(10500)
a1.withdraw(3500)
a2.deposit(5600)
a3.withdraw(60000)
a3.deposit(1000000)

a1.dump()
a2.dump()
a3.dump()
