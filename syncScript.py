#!/usr/bin/python
from subprocess import call
from os.path import expanduser
from os import chdir
from optparse import OptionParser
import paramiko
import cmd
import re
import time

defaultUsername = "cac112"
defaultPkey = "id_rsa"

def declareGlobals():
  global defaultUsername
  global defaultPkey

class ShellRunner(object):

  def __init__(self):
    self.hosts = []
    self.clients = []
    self.sessions = []

  def add_host(self, args):
    """add_host 
    Add the host to the host list"""
    if args:
      self.hosts.append(args.split(','))
    else:
      print "usage: host <host>,<username>,<path rsa>"

  def connect(self):
    """Connect to all hosts in the hosts list"""
    for host in self.hosts:
      client = paramiko.SSHClient()
      client.load_system_host_keys()
      client.set_missing_host_key_policy(
        paramiko.AutoAddPolicy())
      client.get_host_keys()
      try:
        client.connect(host[0], 
          username=host[1], 
          key_filename=host[2])
        self.clients.append(client)
      except:
        print "Failed to connect to %s" % (host)
        break

  def make_sessions(self):
    for (idx, client) in enumerate(self.clients):
      try:
        session = client.get_transport().open_session()
        session.setblocking(blocking=0)
        session.get_pty()
        session.invoke_shell()
        self.sessions.append(session)
        print "Created a session for host %s" % (self.hosts[idx][0])
      except:
        print "Failed to establish session for host %s" % (self.hosts[idx][0])
        break

  def runCommands(self, commands):
    for host, session in zip(self.hosts, self.sessions):
      for command in commands:
        print command
        session.send(command + "\n")
        time.sleep(1)
        stdout=""
        while (session.recv_ready()):
          stdout = stdout + session.recv(2048)
          time.sleep(0.1)
      session.close()

      print "Ran commands for %s" % (host[0])

  def close(self):
    for conn in self.connections:
      conn.close()

def command_line_arguments(parser):
  parser.add_option("--user", dest="username", action="store",
    type="string", default=defaultUsername)

  parser.add_option("--pkey", dest="pkey", action="store",
    type="string", default=defaultPkey)

  parser.add_option("--pass", dest="password", action="store", type="string",
    help="password")

  return parser.parse_args()

def main():
  declareGlobals()
  parser = OptionParser(usage="")
  (options, args) = command_line_arguments(parser)

  username = options.username
  pkey = options.pkey
  if not options.password:
    parser.error("You are required to provide a password field")
  password = options.password

  """sh = ShellRunner()
  with open("hosts", "r") as f:
    for host in f:
      sh.add_host("%s,%s,%s" % (host.rstrip(), username, pkey))
  sh.connect()
  sh.make_sessions()

  commands = ["cd mengIndividual", "git stash", "git pull", username, password]
  commands.extend(["cd code/algorithms", "make compile_multinode"])
  sh.runCommands(commands)"""

  client = paramiko.SSHClient()
  client.load_system_host_keys()
  client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
  client.get_host_keys()
  client.connect("cloud-vm-47-61.doc.ic.ac.uk", username=username, key_filename=pkey)
  
  ftp = client.open_sftp()
  ftp.put('input/10D_500K_synthetic.ds', '/home/cac112/homedir/cac112/input/10D_500K_synthetic.ds')
  ftp.close()

if __name__ == "__main__":
  main()
