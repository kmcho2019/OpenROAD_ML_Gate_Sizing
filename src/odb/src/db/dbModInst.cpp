///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2020, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Generator Code Begin Cpp
#include "dbModInst.h"

#include "dbBlock.h"
#include "dbDatabase.h"
#include "dbDiff.hpp"
#include "dbHashTable.hpp"
#include "dbModITerm.h"
#include "dbModule.h"
#include "dbTable.h"
#include "dbTable.hpp"
#include "odb/db.h"
// User Code Begin Includes
#include "dbGroup.h"
#include "dbModBTerm.h"
#include "dbModuleModInstModITermItr.h"
// User Code End Includes
namespace odb {
template class dbTable<_dbModInst>;

bool _dbModInst::operator==(const _dbModInst& rhs) const
{
  if (_name != rhs._name) {
    return false;
  }
  if (_next_entry != rhs._next_entry) {
    return false;
  }
  if (_parent != rhs._parent) {
    return false;
  }
  if (_module_next != rhs._module_next) {
    return false;
  }
  if (_master != rhs._master) {
    return false;
  }
  if (_group_next != rhs._group_next) {
    return false;
  }
  if (_group != rhs._group) {
    return false;
  }
  if (_moditerms != rhs._moditerms) {
    return false;
  }

  return true;
}

bool _dbModInst::operator<(const _dbModInst& rhs) const
{
  // User Code Begin <
  if (strcmp(_name, rhs._name) >= 0) {
    return false;
  }
  // User Code End <
  return true;
}

void _dbModInst::differences(dbDiff& diff,
                             const char* field,
                             const _dbModInst& rhs) const
{
  DIFF_BEGIN
  DIFF_FIELD(_name);
  DIFF_FIELD(_next_entry);
  DIFF_FIELD(_parent);
  DIFF_FIELD(_module_next);
  DIFF_FIELD(_master);
  DIFF_FIELD(_group_next);
  DIFF_FIELD(_group);
  DIFF_FIELD(_moditerms);
  DIFF_END
}

void _dbModInst::out(dbDiff& diff, char side, const char* field) const
{
  DIFF_OUT_BEGIN
  DIFF_OUT_FIELD(_name);
  DIFF_OUT_FIELD(_next_entry);
  DIFF_OUT_FIELD(_parent);
  DIFF_OUT_FIELD(_module_next);
  DIFF_OUT_FIELD(_master);
  DIFF_OUT_FIELD(_group_next);
  DIFF_OUT_FIELD(_group);
  DIFF_OUT_FIELD(_moditerms);

  DIFF_END
}

_dbModInst::_dbModInst(_dbDatabase* db)
{
  // User Code Begin Constructor
  _name = nullptr;
  _parent = 0;
  _module_next = 0;
  _moditerms = 0;
  _master = 0;
  _group = 0;
  _group_next = 0;
  // User Code End Constructor
}

_dbModInst::_dbModInst(_dbDatabase* db, const _dbModInst& r)
{
  _name = r._name;
  _next_entry = r._next_entry;
  _parent = r._parent;
  _module_next = r._module_next;
  _master = r._master;
  _group_next = r._group_next;
  _group = r._group;
  _moditerms = r._moditerms;
}

dbIStream& operator>>(dbIStream& stream, _dbModInst& obj)
{
  stream >> obj._name;
  stream >> obj._next_entry;
  stream >> obj._parent;
  stream >> obj._module_next;
  stream >> obj._master;
  stream >> obj._group_next;
  stream >> obj._group;
  // User Code Begin >>
  dbBlock* block = (dbBlock*) (obj.getOwner());
  _dbDatabase* db = (_dbDatabase*) (block->getDataBase());
  if (db->isSchema(db_schema_update_hierarchy)) {
    stream >> obj._moditerms;
  }
  // User Code End >>
  return stream;
}

dbOStream& operator<<(dbOStream& stream, const _dbModInst& obj)
{
  stream << obj._name;
  stream << obj._next_entry;
  stream << obj._parent;
  stream << obj._module_next;
  stream << obj._master;
  stream << obj._group_next;
  stream << obj._group;
  // User Code Begin <<
  dbBlock* block = (dbBlock*) (obj.getOwner());
  _dbDatabase* db = (_dbDatabase*) (block->getDataBase());
  if (db->isSchema(db_schema_update_hierarchy)) {
    stream << obj._moditerms;
  }
  // User Code End <<
  return stream;
}

_dbModInst::~_dbModInst()
{
  if (_name) {
    free((void*) _name);
  }
}

////////////////////////////////////////////////////////////////////
//
// dbModInst - Methods
//
////////////////////////////////////////////////////////////////////

const char* dbModInst::getName() const
{
  _dbModInst* obj = (_dbModInst*) this;
  return obj->_name;
}

dbModule* dbModInst::getParent() const
{
  _dbModInst* obj = (_dbModInst*) this;
  if (obj->_parent == 0) {
    return nullptr;
  }
  _dbBlock* par = (_dbBlock*) obj->getOwner();
  return (dbModule*) par->_module_tbl->getPtr(obj->_parent);
}

dbModule* dbModInst::getMaster() const
{
  _dbModInst* obj = (_dbModInst*) this;
  if (obj->_master == 0) {
    return nullptr;
  }
  _dbBlock* par = (_dbBlock*) obj->getOwner();
  return (dbModule*) par->_module_tbl->getPtr(obj->_master);
}

dbGroup* dbModInst::getGroup() const
{
  _dbModInst* obj = (_dbModInst*) this;
  if (obj->_group == 0) {
    return nullptr;
  }
  _dbBlock* par = (_dbBlock*) obj->getOwner();
  return (dbGroup*) par->_group_tbl->getPtr(obj->_group);
}

// User Code Begin dbModInstPublicMethods
dbModInst* dbModInst::create(dbModule* parentModule,
                             dbModule* masterModule,
                             const char* name)
{
  _dbModule* module = (_dbModule*) parentModule;
  _dbBlock* block = (_dbBlock*) module->getOwner();
  _dbModule* master = (_dbModule*) masterModule;

  if (master->_mod_inst != 0) {
    return nullptr;
  }

  dbModInst* ret = nullptr;
  ret = ((dbModule*) module)->findModInst(name);
  if (ret) {
    return nullptr;
  }

  _dbModInst* modinst = block->_modinst_tbl->create();
  modinst->_name = strdup(name);
  ZALLOCATED(modinst->_name);
  modinst->_master = master->getOID();
  modinst->_parent = module->getOID();
  // push to head of list in block
  modinst->_module_next = module->_modinsts;
  module->_modinsts = modinst->getOID();
  master->_mod_inst = modinst->getOID();
  block->_modinst_hash.insert(modinst);
  return (dbModInst*) modinst;
}

void dbModInst::destroy(dbModInst* modinst)
{
  _dbModInst* _modinst = (_dbModInst*) modinst;
  _dbBlock* block = (_dbBlock*) _modinst->getOwner();
  _dbModule* module = (_dbModule*) modinst->getParent();

  _dbModule* master = (_dbModule*) modinst->getMaster();
  master->_mod_inst = dbId<_dbModInst>();  // clear
  dbModule::destroy((dbModule*) master);

  // remove the moditerm connections
  for (auto moditerm : modinst->getModITerms()) {
    moditerm->disconnect();
  }
  // remove the moditerms
  for (auto moditerm : modinst->getModITerms()) {
    block->_moditerm_tbl->destroy((_dbModITerm*) moditerm);
  }
  // unlink from parent start
  uint id = _modinst->getOID();
  _dbModInst* prev = nullptr;
  uint cur = module->_modinsts;
  while (cur) {
    _dbModInst* c = block->_modinst_tbl->getPtr(cur);
    if (cur == id) {
      if (prev == nullptr) {
        module->_modinsts = _modinst->_module_next;
      } else {
        prev->_module_next = _modinst->_module_next;
      }
      break;
    }
    prev = c;
    cur = c->_module_next;
  }
  // unlink from parent end
  if (_modinst->_group) {
    modinst->getGroup()->removeModInst(modinst);
  }
  dbProperty::destroyProperties(_modinst);
  block->_modinst_hash.remove(_modinst);
  block->_modinst_tbl->destroy(_modinst);
}

dbSet<dbModInst>::iterator dbModInst::destroy(dbSet<dbModInst>::iterator& itr)
{
  dbModInst* modinst = *itr;
  dbSet<dbModInst>::iterator next = ++itr;
  destroy(modinst);
  return next;
}

dbSet<dbModITerm> dbModInst::getModITerms()
{
  _dbModInst* _mod_inst = (_dbModInst*) this;
  _dbBlock* _block = (_dbBlock*) _mod_inst->getOwner();
  return dbSet<dbModITerm>(_mod_inst, _block->_module_modinstmoditerm_itr);
}

dbModInst* dbModInst::getModInst(dbBlock* block_, uint dbid_)
{
  _dbBlock* block = (_dbBlock*) block_;
  return (dbModInst*) block->_modinst_tbl->getPtr(dbid_);
}

std::string dbModInst::getHierarchicalName() const
{
  _dbModInst* _obj = (_dbModInst*) this;
  dbBlock* block = (dbBlock*) _obj->getOwner();
  std::string inst_name = std::string(getName());
  dbModule* parent = getParent();
  if (parent == block->getTopModule()) {
    return inst_name;
  }
  return parent->getModInst()->getHierarchicalName() + "/" + inst_name;
}

dbModITerm* dbModInst::findModITerm(const char* name)
{
  dbSet<dbModITerm> moditerms = getModITerms();
  for (dbModITerm* mod_iterm : moditerms) {
    if (!strcmp(mod_iterm->getName(), name)) {
      return mod_iterm;
    }
  }
  return nullptr;
}

void dbModInst::RemoveUnusedPortsAndPins()
{
  std::set<dbModITerm*> kill_set;
  dbModule* module = this->getMaster();
  dbSet<dbModITerm> moditerms = getModITerms();
  dbSet<dbModBTerm> modbterms = module->getModBTerms();

  for (dbModITerm* mod_iterm : moditerms) {
    dbModBTerm* mod_bterm = module->findModBTerm(mod_iterm->getName());
    dbModNet* mod_net = mod_bterm->getModNet();
    dbSet<dbModITerm> dest_mod_iterms = mod_net->getModITerms();
    dbSet<dbBTerm> dest_bterms = mod_net->getBTerms();
    dbSet<dbITerm> dest_iterms = mod_net->getITerms();
    if (dest_mod_iterms.size() == 0 && dest_bterms.size() == 0
        && dest_iterms.size() == 0) {
      kill_set.insert(mod_iterm);
    }
  }
  moditerms = getModITerms();
  modbterms = module->getModBTerms();
  for (auto mod_iterm : kill_set) {
    dbModBTerm* mod_bterm = module->findModBTerm(mod_iterm->getName());
    dbModNet* modbterm_m_net = mod_bterm->getModNet();
    mod_bterm->disconnect();
    dbModBTerm::destroy(mod_bterm);
    dbModNet* moditerm_m_net = mod_iterm->getModNet();
    mod_iterm->disconnect();
    dbModITerm::destroy(mod_iterm);
    if (modbterm_m_net && modbterm_m_net->getBTerms().size() == 0
        && modbterm_m_net->getITerms().size() == 0
        && modbterm_m_net->getModITerms().size() == 0
        && modbterm_m_net->getModBTerms().size() == 0) {
      dbModNet::destroy(modbterm_m_net);
    }
    if (moditerm_m_net && moditerm_m_net->getBTerms().size() == 0
        && moditerm_m_net->getITerms().size() == 0
        && moditerm_m_net->getModITerms().size() == 0
        && moditerm_m_net->getModBTerms().size() == 0) {
      dbModNet::destroy(moditerm_m_net);
    }
  }
}

// User Code End dbModInstPublicMethods
}  // namespace odb
// Generator Code End Cpp
