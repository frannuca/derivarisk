kafka-topics --zookeeper localhost:2181 --topic analytics --create --partitions 3 --replication-factor 1
kafka-topics --zookeeper localhost:2181 --list
kafka-topics --zookeeper localhost:2181 --topic analyticsX --delete