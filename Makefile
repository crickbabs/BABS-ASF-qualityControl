#
# nourdine.bah@crick.ac.uk
#

CC = groovyc
ARCH = jar
OPT = cvf
PACKAGE = qcpipeline

.PHONY: .FORCE
.FORCE:

all: jar

jar: .FORCE $(PACKAGE)/Project.class
	$(ARCH) $(OPT) $(PACKAGE).jar $(PACKAGE)
	rm -v $(PACKAGE)/*.class

$(PACKAGE)/Project.class: \
	.FORCE \
	$(PACKAGE)/Sample.class \
	$(PACKAGE)/SampleSet.class
	$(CC) -d . $(PACKAGE)/Project.gy

$(PACKAGE)/SampleSet.class: .FORCE
	$(CC) -d . $(PACKAGE)/SampleSet.gy

$(PACKAGE)/Sample.class: .FORCE $(PACKAGE)/Fastq.class
	$(CC) -d . $(PACKAGE)/Sample.gy

$(PACKAGE)/Fastq.class: \
	.FORCE \
	$(PACKAGE)/Run.class \
	$(PACKAGE)/Machine.class \
	$(PACKAGE)/FlowCell.class \
	$(PACKAGE)/Lane.class \
	$(PACKAGE)/Read.class
	$(CC) -d . $(PACKAGE)/Fastq.gy

$(PACKAGE)/Run.class: .FORCE
	$(CC) -d . $(PACKAGE)/Run.gy

$(PACKAGE)/Machine.class: .FORCE
	$(CC) -d . $(PACKAGE)/Machine.gy

$(PACKAGE)/FlowCell.class: .FORCE
	$(CC) -d . $(PACKAGE)/FlowCell.gy

$(PACKAGE)/Lane.class: .FORCE
	$(CC) -d . $(PACKAGE)/Lane.gy

$(PACKAGE)/Read.class: .FORCE
	$(CC) -d . $(PACKAGE)/Read.gy




